# 🧬 Genomic Data Visualization Dashboard
## 📌 Project Overview
This project is an interactive Streamlit-based dashboard designed to download, preprocess, analyze, and visualize gene expression data from the NCBI GEO (Gene Expression Omnibus) database. It demonstrates practical understanding of bioinformatics data handling, dimensionality reduction, and exploratory data analysis (EDA) using Python.

The application allows users to:
- Fetch public GEO datasets using accession IDs
- Preprocess and standardize gene expression matrices
- Perform Principal Component Analysis (PCA)
- Visualize expression patterns interactively
- Download processed datasets for downstream analysis

This project highlights skills relevant to bioinformatics analyst / computational biology roles, especially in working with real-world genomic datasets.

---
## 🎯 Key Objectives
- Access and parse large-scale genomic datasets programmatically
- Apply statistical preprocessing techniques to gene expression data
- Reduce dimensionality for biological interpretation
- Create an intuitive and interactive visualization interface
- Ensure reproducibility and usability for researchers

---
## 🛠️ Technologies & Libraries Used
|Category|Tools / Libraries|
|---|---|
|Web App Framework|Streamlit|
|Data Handling|Pandas, NumPy|
|Visualization|Plotly Express, Plotly Graph Objects|
|Machine Learning|scikit-learn (PCA, StandardScaler)|
|Genomic Data Access|GEOparse|
|Deployment Ready|Streamlit Page Config|

---
## 📂 Project Structure
```bash
streamlitapp.py        # Main Streamlit application
README.md              # Project documentation
```
---
## ⚙️ Application Workflow
### 1️⃣ User Input (GEO Accession ID)
- Users provide a valid GEO accession ID (e.g., GSEXXXX).
- The app validates the input before initiating data retrieval.
### 2️⃣ Data Retrieval from GEO
- GEO datasets are downloaded using GEOparse.
- Expression matrices are extracted from GEO Series Matrix files.
- Metadata is parsed for better biological context.
Why this matters:
> Demonstrates ability to work with public genomic repositories used in real research.
### 3️⃣ Data Preprocessing
- Missing values are handled appropriately
- Expression data is transformed into a clean Pandas DataFrame
- Features are standardized using StandardScaler

```python
scaler = StandardScaler()
scaled_data = scaler.fit_transform(expression_df)
```
Purpose:
> Ensures meaningful PCA results and avoids scale-driven bias.
### 4️⃣ Principal Component Analysis (PCA)
- PCA is applied to reduce high-dimensional gene expression data
- Users can visually explore variance across samples
```python
pca = PCA(n_components=2)
pca_result = pca.fit_transform(scaled_data)
```
Biological relevance:
> PCA helps identify clustering patterns, batch effects, or sample similarities.
### 5️⃣ Interactive Visualizations
- Scatter plots for PCA components
- Dynamic hover tooltips for data interpretation
- Responsive layout for better usability
Powered by Plotly, ensuring:
- Zoom
- Pan
- Interactive legends
### 6️⃣ Data Export
- Processed expression data can be downloaded as a CSV file
```python
st.sidebar.download_button(
    label="Download CSV",
    data=csv,
    file_name="expression_data.csv",
    mime="text/csv"
)
```
Use case:
> Enables downstream analysis in R, Python, or other bioinformatics tools.

## 🖥️ User Interface Features
- Wide-layout dashboard for better visualization
- Sidebar controls for inputs and downloads
- Real-time feedback using spinners and status messages
- Clean footer with data source attribution

## 📊 Data Source
- NCBI Gene Expression Omnibus (GEO)
- Public, peer-reviewed genomic datasets

## 🚀 How to Run the Project Locally
### Prerequisites
- Python ≥ 3.8
**Installation**
```bash
pip install streamlit pandas numpy plotly scikit-learn GEOparse
```
**Run the App**
```bash
streamlit run streamlitapp.py
```
## 🧪 Example Use Cases
- Exploratory analysis of RNA-seq / microarray datasets
- Teaching PCA and dimensionality reduction in genomics
- Rapid hypothesis generation before differential expression analysis

## 🔮 Future Improvements
- Differential gene expression analysis
- Sample annotation & phenotype integration
- Additional dimensionality reduction methods (t-SNE, UMAP)
- Cloud deployment (Streamlit Cloud)

## 👩‍🔬 Author
Shelinpreethi
Aspiring Bioinformatics Analyst

## 📎 Acknowledgements
NCBI GEO Database
Streamlit Open Source Community
