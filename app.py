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

        # Constructing the query for PubMed search
        if start_year and end_year:
            query_with_dates = f"{user_query} AND ({start_year}[PD] : {end_year}[PD])"
        elif start_year:
            query_with_dates = f"{user_query} AND ({start_year}[PD])"
        elif end_year:
            query_with_dates = f"{user_query} AND ({end_year}[PD])"
        else:
            query_with_dates = user_query

        results_pubmed = search_pubmed(query_with_dates)
        results_scholarly = search_scholarly(user_query)

    return render_template('index.html', user=user, results_pubmed=results_pubmed, results_scholarly=results_scholarly)

@app.route('/register', methods=['GET', 'POST'])
def register():
    if request.method == 'POST':
        email = request.form['email']
        password = request.form['password']
        try:
            user = auth.create_user(email=email, password=password)
            flash("Inscription réussie. Vous pouvez maintenant vous connecter.")
            return redirect(url_for('login'))
        except Exception as e:
            flash(f"Erreur d'inscription : {e}")
    return render_template('register.html')

@app.route('/login', methods=['GET', 'POST'])
def login():
    if request.method == 'POST':
        email = request.form['email']
        password = request.form['password']
        try:
            user = auth.get_user_by_email(email)
            # Use Firebase Admin SDK to verify the password (this requires additional setup)
            # Unfortunately, Firebase Admin SDK does not provide a way to sign in with email and password.
            # You'll need to use Firebase Client SDK or implement your own verification.
            session["user"] = {
                "id": user.uid,
                "email": user.email
            }
            flash("Connexion réussie.")
            return redirect(url_for('index'))
        except Exception as e:
            flash(f"Erreur de connexion : {e}")
    return render_template('login.html')

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
