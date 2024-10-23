import firebase_admin
from firebase_admin import credentials, auth, firestore
from flask import Flask, render_template, request, redirect, url_for, session, flash
import pandas as pd
from Bio import Entrez
from scholarly import scholarly
import requests

app = Flask(__name__)
app.secret_key = 'your_secret_key'  # Remplacez par une clé secrète sécurisée.

Entrez.email = "your_email@example.com"  # Remplacez par votre email.

# Firebase & Firestore initialization
cred = credentials.Certificate("serviceAccountKey.json")
firebase_admin.initialize_app(cred)
db = firestore.client()

# Identifiants OAuth (Firebase Console)
GOOGLE_CLIENT_ID = "901804151640-nqv187k3fvfq9d4rnvf5dftsd2fqhjvb.apps.googleusercontent.com"
GOOGLE_CLIENT_SECRET = "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"

def handle_pubmed_search(query):
    """Recherche PubMed avec gestion des erreurs."""
    try:
        handle = Entrez.esearch(db="pubmed", term=query)
        record = Entrez.read(handle)
        handle.close()

        ids = record.get("IdList", [])
        if not ids:
            return []

        id_str = ",".join(ids)
        with Entrez.efetch(db="pubmed", id=id_str, retmode="xml") as handle:
            articles = Entrez.read(handle)

        return [
            {
                "Titre": article.get("MedlineCitation", {}).get("Article", {}).get("ArticleTitle", "Titre non disponible"),
                "Année": article.get("MedlineCitation", {}).get("Article", {}).get("ArticleDate", [{}])[0].get("Year", "N/A"),
                "Lien": f"https://pubmed.ncbi.nlm.nih.gov/{article['MedlineCitation'].get('PMID', 'N/A')}/"
            }
            for article in articles.get("PubmedArticle", [])
        ]
    except Exception as e:
        print(f"Erreur lors de la recherche PubMed : {e}")
        return []

def handle_scholarly_search(query):
    """Recherche avec Scholarly."""
    results = []
    try:
        search_query = scholarly.search_pubs(query)
        for _ in range(10):
            pub = next(search_query)
            results.append({
                "Titre": pub.get("bib", {}).get("title", "Titre non disponible"),
                "Année": pub.get("bib", {}).get("pub_year", "Année non disponible"),
                "Lien": pub.get("pub_url", "Non disponible")
            })
    except (StopIteration, Exception) as e:
        print(f"Erreur lors de la recherche Scholarly : {e}")
    return results

def ajouter_article_favori(user_id, article):
    """Ajout d'un article aux favoris de l'utilisateur."""
    db.collection('users').document(user_id).collection('panier').document(article['Titre']).set(article)

@app.route('/', methods=['GET', 'POST'])
def index():
    user = session.get("user")
    results_pubmed, results_scholarly = [], []

    if request.method == 'POST':
        user_query = request.form.get('query')
        start_year, end_year = request.form.get('start_year'), request.form.get('end_year')

        query = f"{user_query} AND ({start_year}[PD] : {end_year}[PD])" if start_year and end_year else user_query
        results_pubmed = handle_pubmed_search(query)
        results_scholarly = handle_scholarly_search(user_query)

    return render_template('index.html', user=user, results_pubmed=results_pubmed, results_scholarly=results_scholarly)

@app.route('/register', methods=['GET', 'POST'])
def register():
    if request.method == 'POST':
        email, password = request.form['email'], request.form['password']
        try:
            auth.create_user(email=email, password=password)
            flash("Utilisateur enregistré avec succès.")
            return redirect(url_for('login'))
        except Exception as e:
            flash(f"Erreur d'enregistrement : {e}")
    return render_template('register.html')

@app.route('/login', methods=['GET', 'POST'])
def login():
    if request.method == 'POST':
        email, password = request.form['email'], request.form['password']
        try:
            user = auth.get_user_by_email(email)
            session['user'] = {"id": user.uid, "email": email}
            return redirect(url_for('index'))
        except Exception as e:
            flash(f"Erreur de connexion : {e}")
    return render_template('login.html')

@app.route('/logout')
def logout():
    session.pop('user', None)
    flash("Vous avez été déconnecté.")
    return redirect(url_for('index'))

@app.route('/ajouter_article', methods=['POST'])
def ajouter_article():
    if "user" not in session:
        flash("Veuillez vous connecter.")
        return redirect(url_for("login"))

    article = {key: request.form[key] for key in ["titre", "annee", "lien"]}
    ajouter_article_favori(session["user"]["id"], article)
    flash("Article ajouté aux favoris.")
    return redirect(url_for("index"))

@app.route('/panier')
def panier():
    user = session.get("user")
    if not user:
        flash("Veuillez vous connecter.")
        return redirect(url_for("login"))

    articles = [
        doc.to_dict() for doc in db.collection('users').document(user["id"]).collection('panier').stream()
    ]
    return render_template('panier.html', user=user, articles=articles)

if __name__ == '__main__':
    app.run(debug=True)
