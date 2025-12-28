import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import numpy as np
import os
import GEOparse
import time

# Page configuration
st.set_page_config(
    page_title="Genomic Data Visualization Dashboard",
    page_icon="🧬",
    layout="wide"
)
linkedin_url = "https://www.linkedin.com/in/shelinpreethi/"
github_url = "genomic-data-visualization-dashboard"

# Title and description
st.title("🧬 Genomic Data Visualization Dashboard")
st.markdown("Developed by [Preethi Shelin]({})".format(linkedin_url))
st.markdown("""
This dashboard visualizes gene expression patterns from public genomic datasets.
Explore gene expression across samples, perform PCA analysis, and identify differentially expressed genes.
""")

# Sidebar for user inputs
st.sidebar.header("GEO Data Input")
geo_accession_id = st.sidebar.text_input("Enter GEO Accession ID (e.g., GSE5364, GSE48350, GSE10072)", "GSE5364")


@st.cache_data
def load_data(accession_id):
    """
    Load genomic data from GEO using GEOparse.
    Properly extracts expression data from GSM objects within the GSE.
    """
    gse = None
    local_filepath = os.path.join("./GEO_data", f"{accession_id}_family.soft.gz")

    if not accession_id:
        st.warning("Please provide a GEO Accession ID.")
        return pd.DataFrame(), pd.DataFrame()

    # Attempt to download from GEO
    status_placeholder = st.empty()
    status_placeholder.info(f"Attempting to download/load data for {accession_id} from GEO...")
    
    try:
        gse = GEOparse.get_GEO(geo=accession_id, destdir="./GEO_data", silent=False)
        status_placeholder.success(f"Successfully downloaded and parsed {accession_id} from GEO.")
        time.sleep(1)
        status_placeholder.empty()
    except Exception as e:
        status_placeholder.warning(f"Automated download from GEO failed: {e}. Trying local file...")
        time.sleep(1)
        
        # Fallback to local file
        if os.path.exists(local_filepath):
            try:
                gse = GEOparse.get_GEO(filepath=local_filepath, silent=False)
                status_placeholder.success(f"Successfully loaded local file: {local_filepath}")
                time.sleep(1)
                status_placeholder.empty()
            except Exception as parse_e:
                status_placeholder.error(f"Failed to parse local file '{local_filepath}': {parse_e}")
                return pd.DataFrame(), pd.DataFrame()
        else:
            status_placeholder.error(f"Download failed AND local file '{local_filepath}' not found.")
            st.error("Please check your network connection or manually download the SOFT file.")
            return pd.DataFrame(), pd.DataFrame()

    if not gse:
        return pd.DataFrame(), pd.DataFrame()
    
    # --- Extract expression data from GSM objects ---
    try:
        status_placeholder.info(f"Found {len(gse.gsms)} samples in the dataset. Processing...")
        
        # Collect expression data from all samples
        expression_dict = {}
        sample_ids = []
        
        for gsm_name, gsm in gse.gsms.items():
            # Check if the sample has a table with expression data
            if hasattr(gsm, 'table') and not gsm.table.empty:
                # Get the gene identifier column
                gene_id_col = None
                possible_cols = ['ID_REF', 'ID', 'PROBE_ID', 'ProbeID']
                for col in possible_cols:
                    if col in gsm.table.columns:
                        gene_id_col = col
                        break
                
                # Get the value column
                value_col = None
                possible_value_cols = ['VALUE', 'value', 'Value']
                for col in possible_value_cols:
                    if col in gsm.table.columns:
                        value_col = col
                        break
                
                if gene_id_col and value_col:
                    # Store expression values for this sample
                    sample_data = gsm.table.set_index(gene_id_col)[value_col]
                    expression_dict[gsm_name] = sample_data
                    sample_ids.append(gsm_name)
        
        if not expression_dict:
            status_placeholder.error("No valid expression data found in any samples.")
            return pd.DataFrame(), pd.DataFrame()
        
        # Create expression dataframe
        expression_df = pd.DataFrame(expression_dict)
        expression_df.index.name = 'Gene'
        
        # Convert to numeric, replacing non-numeric values with NaN
        expression_df = expression_df.apply(pd.to_numeric, errors='coerce')
        
        # Remove rows with too many missing values
        expression_df = expression_df.dropna(thresh=len(expression_df.columns) * 0.5)
        
        status_placeholder.success(f"Expression matrix created: {expression_df.shape[0]} genes × {expression_df.shape[1]} samples")
        time.sleep(1)
        status_placeholder.empty()
        
    except Exception as e:
        status_placeholder.error(f"Error extracting expression data: {e}")
        import traceback
        st.error(traceback.format_exc())
        return pd.DataFrame(), pd.DataFrame()
    
    # --- Extract metadata ---
    try:
        metadata_list = []
        for gsm_name, gsm in gse.gsms.items():
            if gsm_name in expression_df.columns:
                metadata_dict = {
                    'Sample': gsm_name,
                    'title': gsm.metadata.get('title', [''])[0] if 'title' in gsm.metadata else '',
                    'source': gsm.metadata.get('source_name_ch1', [''])[0] if 'source_name_ch1' in gsm.metadata else '',
                }
                
                # Try to extract characteristics
                if 'characteristics_ch1' in gsm.metadata:
                    for char in gsm.metadata['characteristics_ch1']:
                        if ':' in char:
                            key, value = char.split(':', 1)
                            metadata_dict[key.strip()] = value.strip()
                
                metadata_list.append(metadata_dict)
        
        metadata_df = pd.DataFrame(metadata_list)
        metadata_df = metadata_df.set_index('Sample')
        
        # Try to identify condition column
        condition_col = None
        possible_condition_cols = ['condition', 'group', 'treatment', 'disease', 'tissue']
        
        for col in metadata_df.columns:
            col_lower = col.lower()
            if any(keyword in col_lower for keyword in possible_condition_cols):
                condition_col = col
                break
        
        if condition_col:
            metadata_df['Condition'] = metadata_df[condition_col]
        elif 'source' in metadata_df.columns:
            metadata_df['Condition'] = metadata_df['source']
        else:
            # Create a default condition based on sample names
            metadata_df['Condition'] = 'Sample'
        
        # Reset index to have Sample as a column too
        metadata_df = metadata_df.reset_index()
        
    except Exception as e:
        status_placeholder.error(f"Error extracting metadata: {e}")
        import traceback
        st.error(traceback.format_exc())
        return pd.DataFrame(), pd.DataFrame()
    
    return expression_df, metadata_df


# Load data based on user input
if geo_accession_id:
    expression_df, metadata_df = load_data(geo_accession_id)
else:
    expression_df, metadata_df = pd.DataFrame(), pd.DataFrame()

# Only proceed if data is loaded successfully
if not expression_df.empty and not metadata_df.empty:
    st.sidebar.markdown("---")
    analysis_type = st.sidebar.selectbox(
        "Select Analysis",
        ["Gene Expression Heatmap", "PCA Analysis", "Volcano Plot", "Box Plot Comparison"]
    )

    # Data preview section
    with st.expander("📊 View Raw Data"):
        st.subheader("Expression Data Preview")
        st.dataframe(expression_df.head(10))
        st.subheader("Sample Metadata")
        st.dataframe(metadata_df)

    # Main analysis section
    if analysis_type == "Gene Expression Heatmap":
        st.header("Gene Expression Heatmap")
        
        col1, col2 = st.columns([1, 3])
        
        with col1:
            max_genes = min(100, len(expression_df))
            n_genes_display = st.slider("Number of genes to display", 10, max_genes, min(50, max_genes))
            clustering = st.checkbox("Apply clustering", value=True)
        
        with col2:
            gene_variance = expression_df.var(axis=1).sort_values(ascending=False)
            top_genes = gene_variance.head(n_genes_display).index
            
            heatmap_data = expression_df.loc[top_genes]
            
            fig = px.imshow(
                heatmap_data,
                labels=dict(x="Samples", y="Genes", color="Expression"),
                aspect="auto",
                color_continuous_scale="RdBu_r"
            )
            fig.update_layout(height=600)
            st.plotly_chart(fig, use_container_width=True)

    elif analysis_type == "PCA Analysis":
        st.header("Principal Component Analysis (PCA)")
        
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(expression_df.T)
        
        n_components = min(10, scaled_data.shape[0], scaled_data.shape[1])
        pca = PCA(n_components=n_components)
        pca_result = pca.fit_transform(scaled_data)
        
        pca_df = pd.DataFrame(
            pca_result[:, :min(3, n_components)],
            columns=[f'PC{i+1}' for i in range(min(3, n_components))],
            index=expression_df.columns
        )
        
        # Keep Sample as a column by resetting index before join
        pca_df['Sample'] = pca_df.index
        pca_df = pca_df.merge(metadata_df, on='Sample', how='left')
        
        col1, col2 = st.columns(2)
        
        with col1:
            fig = px.scatter(
                pca_df,
                x='PC1',
                y='PC2',
                color='Condition',
                hover_data=['Sample'],
                title=f"PCA: PC1 vs PC2<br>Variance Explained: PC1={pca.explained_variance_ratio_[0]:.2%}, PC2={pca.explained_variance_ratio_[1]:.2%}"
            )
            st.plotly_chart(fig, use_container_width=True)
        
        with col2:
            variance_df = pd.DataFrame({
                'PC': [f'PC{i+1}' for i in range(len(pca.explained_variance_ratio_))],
                'Variance Explained': pca.explained_variance_ratio_ * 100
            })
            
            fig = px.bar(
                variance_df,
                x='PC',
                y='Variance Explained',
                title="Scree Plot - Variance Explained by Each PC"
            )
            st.plotly_chart(fig, use_container_width=True)

    elif analysis_type == "Volcano Plot":
        st.header("Volcano Plot - Differential Expression")
        
        # Get unique conditions
        unique_conditions = metadata_df['Condition'].unique()
        
        if len(unique_conditions) >= 2:
            col1, col2 = st.columns(2)
            with col1:
                control_condition = st.selectbox("Select Control Condition", unique_conditions, index=0)
            with col2:
                treatment_condition = st.selectbox("Select Treatment Condition", unique_conditions, index=1 if len(unique_conditions) > 1 else 0)
            
            control_samples = metadata_df[metadata_df['Condition'] == control_condition]['Sample'].values
            treatment_samples = metadata_df[metadata_df['Condition'] == treatment_condition]['Sample'].values
            
            if len(control_samples) > 0 and len(treatment_samples) > 0:
                control_mean = expression_df[control_samples].mean(axis=1)
                treatment_mean = expression_df[treatment_samples].mean(axis=1)
                
                log2_fc = treatment_mean - control_mean
                
                # Mock p-values (in real analysis, use scipy.stats.ttest_ind)
                from scipy import stats
                p_values = []
                for gene in expression_df.index:
                    try:
                        _, p_val = stats.ttest_ind(
                            expression_df.loc[gene, control_samples],
                            expression_df.loc[gene, treatment_samples]
                        )
                        p_values.append(p_val if not np.isnan(p_val) else 1.0)
                    except:
                        p_values.append(1.0)
                
                p_values = np.array(p_values)
                neg_log10_p = -np.log10(p_values + 1e-10)  # Add small value to avoid log(0)
                
                volcano_df = pd.DataFrame({
                    'Gene': expression_df.index,
                    'Log2FC': log2_fc,
                    'NegLog10P': neg_log10_p
                })
                
                fc_threshold = st.sidebar.slider("Log2 Fold Change Threshold", 0.5, 3.0, 1.0)
                p_threshold = st.sidebar.slider("-Log10(p-value) Threshold", 1.0, 5.0, 2.0)
                
                volcano_df['Regulation'] = 'Not Significant'
                volcano_df.loc[(volcano_df['Log2FC'] > fc_threshold) & (volcano_df['NegLog10P'] > p_threshold), 'Regulation'] = 'Up-regulated'
                volcano_df.loc[(volcano_df['Log2FC'] < -fc_threshold) & (volcano_df['NegLog10P'] > p_threshold), 'Regulation'] = 'Down-regulated'
                
                fig = px.scatter(
                    volcano_df,
                    x='Log2FC',
                    y='NegLog10P',
                    color='Regulation',
                    hover_data=['Gene'],
                    title=f"Volcano Plot: {treatment_condition} vs {control_condition}",
                    color_discrete_map={
                        'Up-regulated': 'red',
                        'Down-regulated': 'blue',
                        'Not Significant': 'gray'
                    }
                )
                
                fig.add_hline(y=p_threshold, line_dash="dash", line_color="black", opacity=0.5)
                fig.add_vline(x=fc_threshold, line_dash="dash", line_color="black", opacity=0.5)
                fig.add_vline(x=-fc_threshold, line_dash="dash", line_color="black", opacity=0.5)
                
                st.plotly_chart(fig, use_container_width=True)
                
                col1, col2, col3 = st.columns(3)
                col1.metric("Up-regulated Genes", len(volcano_df[volcano_df['Regulation'] == 'Up-regulated']))
                col2.metric("Down-regulated Genes", len(volcano_df[volcano_df['Regulation'] == 'Down-regulated']))
                col3.metric("Total Significant", len(volcano_df[volcano_df['Regulation'] != 'Not Significant']))
            else:
                st.warning("Need at least one sample in each condition for comparison.")
        else:
            st.warning("Need at least 2 different conditions for volcano plot analysis.")

    else:  # Box Plot Comparison
        st.header("Gene Expression Comparison")
        
        max_default = min(5, len(expression_df.index))
        selected_genes = st.multiselect(
            "Select genes to compare",
            options=expression_df.index.tolist(),
            default=expression_df.index[:max_default].tolist()
        )
        
        if selected_genes:
            plot_data = []
            for gene in selected_genes:
                for sample in expression_df.columns:
                    condition = metadata_df[metadata_df['Sample'] == sample]['Condition'].values
                    if len(condition) > 0:
                        plot_data.append({
                            'Gene': gene,
                            'Expression': expression_df.loc[gene, sample],
                            'Condition': condition[0],
                            'Sample': sample
                        })
            
            plot_df = pd.DataFrame(plot_data)
            
            fig = px.box(
                plot_df,
                x='Gene',
                y='Expression',
                color='Condition',
                title="Gene Expression Comparison by Condition",
                points="all"
            )
            fig.update_layout(height=500)
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("Please select at least one gene to visualize")

    # Download section
    st.sidebar.markdown("---")
    st.sidebar.subheader("📥 Download Data")

    if st.sidebar.button("Download Expression Data"):
        csv = expression_df.to_csv()
        st.sidebar.download_button(
            label="Download CSV",
            data=csv,
            file_name=f"{geo_accession_id}_expression_data.csv",
            mime="text/csv"
        )

# Footer
st.markdown("---")
st.markdown("""
**Data Source:** GEO Database (NCBI)  
**Built with:** Streamlit, Plotly, scikit-learn, GEOparse
""")
st.markdown("**GitHub:** [genomic-data-visualization-dashboard]({})".format(github_url))