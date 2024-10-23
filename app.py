import firebase_admin
from firebase_admin import credentials, auth, firestore
from flask import Flask, render_template, request, redirect, url_for, session, flash
import pandas as pd
from Bio import Entrez
from scholarly import scholarly
import requests

app = Flask(__name__)
app.secret_key = 'your_secret_key'  # Replace with your secret key

# Configure the email for Entrez
Entrez.email = "your_email@example.com"  # Replace with your email

# Initialize Firebase
cred = credentials.Certificate("serviceAccountKey.json")  # Path to your Firebase service account key
firebase_admin.initialize_app(cred)

# Initialize Firestore
db = firestore.client()

# OAuth client ID and secret (replace with your values from Firebase Console)
GOOGLE_CLIENT_ID = "901804151640-nqv187k3fvfq9d4rnvf5dftsd2fqhjvb.apps.googleusercontent.com"
GOOGLE_CLIENT_SECRET = "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"

# Function for PubMed search
def search_pubmed(query):
    try:
        handle = Entrez.esearch(db="pubmed", term=query)
        record = Entrez.read(handle)
        handle.close()
        ids = record["IdList"]
        results = []

        if ids:
            id_str = ','.join(ids)
            handle = Entrez.efetch(db="pubmed", id=id_str, retmode="xml")
            articles = Entrez.read(handle)
            handle.close()

            for article in articles.get("PubmedArticle", []):
                citation = article.get("MedlineCitation", {})
                if citation:
                    title = citation.get("Article", {}).get("ArticleTitle", "Titre non disponible")
                    article_date = citation.get("Article", {}).get("ArticleDate", [])
                    year = article_date[0].get("Year", "Année non disponible") if article_date else "Année non disponible"
                    pubmed_id = citation.get('PMID', "ID non disponible")
                    url = f"https://pubmed.ncbi.nlm.nih.gov/{pubmed_id}/"
                    results.append({"Titre": title, "Année": year, "Lien": url})

        return results
    except Exception as e:
        print(f"Error searching PubMed: {e}")
        return []

# Function for scholarly search
def search_scholarly(query):
    results = []
    search_query = scholarly.search_pubs(query)

    for i in range(10):
        try:
            pub = next(search_query)
            title = pub.get("bib", {}).get("title", "Titre non disponible")
            year = pub.get("bib", {}).get("pub_year", "Année non disponible")
            url = pub.get("pub_url", "Non disponible")
            results.append({"Titre": title, "Année": year, "Lien": url})
        except StopIteration:
            break
        except Exception as e:
            print(f"Error searching Scholarly: {e}")
            break

    return results

# Function to add an article to favorites
def ajouter_article_favori(user_id, article):
    user_ref = db.collection('users').document(user_id)
    panier_ref = user_ref.collection('panier').document(article['Titre'])
    panier_ref.set(article)

@app.route('/', methods=['GET', 'POST'])
def index():
    user = session.get("user")
    results_pubmed = []
    results_scholarly = []

    if request.method == 'POST':
        user_query = request.form.get('query')
        start_year = request.form.get('start_year')
        end_year = request.form.get('end_year')

        # PubMed search
        if start_year and end_year:
            query_with_dates = f"{user_query} AND ({start_year}[PD] : {end_year}[PD])"
        elif start_year:
            query_with_dates = f"{user_query} AND ({start_year}[PD])"
        elif end_year:
            query_with_dates = f"{user_query} AND ({end_year}[PD])"
        else:
            query_with_dates = user_query

        results_pubmed = search_pubmed(query_with_dates)

        # Scholarly search
        results_scholarly = search_scholarly(user_query)

    return render_template('index.html', user=user, results_pubmed=results_pubmed, results_scholarly=results_scholarly)

@app.route('/login')
def login():
    google_auth_url = (
        f"https://accounts.google.com/o/oauth2/v2/auth?"
        f"client_id={GOOGLE_CLIENT_ID}&"
        f"redirect_uri=https://web-production-35f5e.up.railway.app&"
        f"response_type=code&"
        f"scope=email profile"
    )
    return redirect(google_auth_url)
@app.route('/callback')
def callback():
    code = request.args.get("code")
    token_url = "https://oauth2.googleapis.com/token"
    data = {
        "code": code,
        "client_id": GOOGLE_CLIENT_ID,
        "client_secret": GOOGLE_CLIENT_SECRET,
        "redirect_uri": "https://web-production-35f5e.up.railway.app/callback",
        "grant_type": "authorization_code",
    }
    
    try:
        token_response = requests.post(token_url, data=data).json()
        access_token = token_response.get("access_token")

        user_info = requests.get(
            "https://www.googleapis.com/oauth2/v1/userinfo?alt=json",
            headers={"Authorization": f"Bearer {access_token}"}
        ).json()

        # Save user info in session including profile picture URL
        session["user"] = {
            "id": user_info["id"],
            "name": user_info["name"],
            "email": user_info["email"],
            "picture": user_info["picture"]  # Add profile picture URL
        }

        # Check if user exists in Firestore; if not, create them
        user_ref = db.collection('users').document(user_info['id'])
        if not user_ref.get().exists:
            user_ref.set({"email": user_info["email"], "name": user_info["name"], "picture": user_info["picture"]})

        return redirect(url_for("index"))
    except Exception as e:
        print(f"Error in callback: {e}")
        flash("Une erreur s'est produite lors de l'authentification.")
        return redirect(url_for("index"))

@app.route('/logout')
def logout():
    session.pop("user", None)
    flash("Vous avez été déconnecté.")
    return redirect(url_for("index"))

@app.route('/ajouter_article', methods=['POST'])
def ajouter_article():
    if "user" not in session:
        flash("Veuillez vous connecter pour ajouter un article à votre panier.")
        return redirect(url_for("login"))

    article = {
        "Titre": request.form["titre"],
        "Année": request.form["annee"],
        "Lien": request.form["lien"]
    }
    user_id = session["user"]["id"]
    ajouter_article_favori(user_id, article)
    flash("Article ajouté à votre panier avec succès.")
    return redirect(url_for("index"))
@app.route('/panier', methods=['GET', 'POST'])
def panier():
    user = session.get("user")
    if not user:
        flash("Veuillez vous connecter pour accéder à votre panier.")
        return redirect(url_for("login"))

    user_id = user["id"]
    # Fetch the user's articles from Firestore
    panier_ref = db.collection('users').document(user_id).collection('panier')
    articles = panier_ref.stream()
    articles_list = [article.to_dict() for article in articles]

    return render_template('panier.html', user=user, articles=articles_list)

@app.route('/add_article', methods=['POST'])
def add_article():
    user = session.get("user")
    if not user:
        flash("Veuillez vous connecter pour ajouter des articles.")
        return redirect(url_for("login"))

    user_id = user["id"]
    titre = request.form.get('titre')
    annee = request.form.get('annee')
    lien = request.form.get('lien')

    # Add article to Firestore
    panier_ref = db.collection('users').document(user_id).collection('panier')
    panier_ref.add({
        'Titre': titre,
        'Année': annee,
        'Lien': lien
    })

    flash("Article ajouté avec succès au panier!")
    return redirect(url_for('panier'))

@app.route('/remove_article/<string:titre>', methods=['POST'])
def remove_article(titre):
    user = session.get("user")
    if not user:
        flash("Veuillez vous connecter pour retirer des articles.")
        return redirect(url_for("login"))

    user_id = user["id"]
    panier_ref = db.collection('users').document(user_id).collection('panier')

    # Query to find the article by title
    query = panier_ref.where('Titre', '==', titre).limit(1).stream()
    for article in query:
        panier_ref.document(article.id).delete()

    flash("Article retiré avec succès du panier!")
    return redirect(url_for('panier'))


if __name__ == '__main__':
    app.run(debug=True)