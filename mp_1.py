import streamlit as st
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import matplotlib.pyplot as plt
import io
import py3Dmol
import requests
import streamlit.components.v1 as components
import base64


# Title
st.title("Protein Sequence Analyzer")


# Sidebar Input
st.sidebar.header("Input")
sequence = st.sidebar.text_area("Enter Protein Sequence")
analyze = st.sidebar.button("Analyze")


if analyze and sequence:
    analysis = ProteinAnalysis(sequence)
    
    st.subheader("Basic Properties")

    col1, col2 = st.columns(2)

    with col1:
        st.write("Molecular Weight:", round(analysis.molecular_weight(), 2))
        st.write("Isoelectric Point:", round(analysis.isoelectric_point(), 2))

    with col2:
        st.write("Aromaticity:", round(analysis.aromaticity(), 2))
        st.write("Instability Index:", round(analysis.instability_index(), 2))

    
    # Amino Acid Composition
  
    st.subheader("Amino Acid Composition")

    aa_percent = analysis.amino_acids_percent

    df = pd.DataFrame({"Amino Acid": list(aa_percent.keys()),"Percentage": list(aa_percent.values())})

    st.dataframe(df)

    
    # Plot
    
    st.subheader("Amino Acid Composition Plot")
    fig1, ax1 = plt.subplots()
    plt.xlabel('X-axis',size=12) 
    plt.ylabel('Y-axis',size=12)
    ax1.bar(df["Amino Acid"], df["Percentage"])
    st.pyplot(fig1)
    
    buf1 = io.BytesIO()
    fig1.savefig(buf1, format="png")
    st.download_button(
        label="Download Composition Plot",
        data=buf1.getvalue(),
        file_name="amino_acid_composition.png",
        mime="image/png"
    )
    
    # Download
    
    st.subheader("Download Results")

    csv = df.to_csv(index=False).encode('utf-8')

    st.download_button(
        label="Download Composition CSV",
        data=csv,
        file_name="protein_analysis.csv",
        mime="text/csv"
    )
    
    
    # Acidic & Basic
    st.subheader("Acidic vs Basic Amino Acids")
    
    D=0
    E=0
    K=0
    R=0
    H=0
    for i in range(0,len(sequence)):
        if(sequence[i]=="D"):
            D=D+1
        elif(sequence[i]=="E"):
            E=E+1
        elif(sequence[i]=="K"):
            K=K+1
        elif(sequence[i]=="R"):
            R=R+1
        elif(sequence[i]=="H"):
            H=H+1
    acidic = ((D+E)/len(sequence))*100
    basic = ((K+R+H)/len(sequence))*100

    st.write("Acidic %:", round(acidic, 3))
    st.write("Basic %:", round(basic, 3))

    fig2, ax2 = plt.subplots()
    ax2.bar(["Acidic", "Basic"], [acidic, basic], color=['red', 'blue'])
    plt.xlabel('X-axis', size=12) 
    plt.ylabel('Y-axis', size=12)  
    ax2.set_title("Acidic vs Basic Content")
    st.pyplot(fig2)
        
    buf2 = io.BytesIO()
    fig2.savefig(buf2, format="png")
    st.download_button(
        label="Download Acidic vs Basic Plot",
        data=buf2.getvalue(),
        file_name="acidic_basic_plot.png",
        mime="image/png"
    )
    
    
    # PTM Prediction
    def predict_ptm_sites(seq):
        ptm_sites = []
        for i, aa in enumerate(seq):
            if aa in ['S', 'T', 'Y']:
                ptm_sites.append((i+1, aa))
        return ptm_sites
    
    st.subheader("Predicted PTM Sites (Phosphorylation)")
    ptm_sites = predict_ptm_sites(sequence)

    if ptm_sites:
        ptm_df = pd.DataFrame(ptm_sites, columns=["Position", "Residue"])
        st.dataframe(ptm_df)
    else:
        st.write("No PTM sites predicted")


st.subheader("Subsequence Finder")

subseq = st.text_input("Enter subsequence to search").upper()

if subseq and sequence:

    positions = []
    for i in range(len(sequence)):
        if sequence[i:i+len(subseq)] == subseq:
            positions.append(i + 1)

    if positions:
        st.success(f"Subsequence found at positions: {positions}")
        st.write(f"Total occurrences: {len(positions)}")
    else:
        st.error("Subsequence not found!")
        
    highlighted = sequence.replace(subseq, f"[{subseq}]")
    st.write("Highlighted Sequence:")
    st.code(highlighted)
    
    
    
st.subheader("3D Protein Structure Viewer")

pdb_id_input = st.text_input("Enter PDB ID for visualization")

if pdb_id_input:
    pdb_id_input = pdb_id_input.upper()


    pdb_url = f"https://files.rcsb.org/download/{pdb_id_input}.pdb"
    pdb_data = requests.get(pdb_url).text

    if "HEADER" in pdb_data:

        st.success(f"PDB {pdb_id_input} loaded!")

        view = py3Dmol.view(width=700, height=500)
        view.addModel(pdb_data, 'pdb')
        view.setStyle({'cartoon': {'color': 'yellow'}})
        view.setBackgroundColor('white')
        view.zoomTo()
        
        components.html(view._make_html(), height=500, width=700)

    else:
        st.error("Invalid PDB ID or structure not found.")

    html_data = view._make_html()
    st.download_button(
        label="Download 3D Structure (HTML)",
        data=html_data,
        file_name=f"{pdb_id_input}_structure.html",
        mime="text/html"
    )
