# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 

# * * * * * * * * FUNCIONES * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  

import streamlit as st
# ~ ~ ~ ~ ~ ~ ~ ~ LIBRERIAS ~ ~ ~ ~ ~ ~ ~ ~

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import os
import json
import zipfile
import shutil
from Bio.PDB import MMCIFParser, PDBIO, Select
from io import BytesIO
import py3Dmol

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

# * * * * * * * * FUNCIONES * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  

# FUNCION- FLOR -----------------------------------------------------------------------------------------------
def load_af3_pdb(zip_file_path):
    name = os.path.basename(zip_file_path).replace('.zip', '')
    folder = f'AF3_files/{name}'
    os.makedirs(folder, exist_ok=True)

    with zipfile.ZipFile(zip_file_path, 'r') as fz:
        fz.extractall(folder)

    pdb_file_path = None
    for path in os.listdir(folder):
        long_path = os.path.join(folder, path)
        if long_path.endswith("_model_0.cif"):
            file_parser = MMCIFParser(QUIET=True)
            structure = file_parser.get_structure("base", long_path)
            io = PDBIO()
            io.set_structure(structure)
            io.save(f'{folder}/{name}_relaxed.pdb', select=Select())
            pdb_file_path = f'{folder}/{name}_relaxed.pdb'
            print(f"Archivo PDB guardado en: {pdb_file_path}")

    return pdb_file_path

def extract_ranking_scores_from_zip(zip_file):
    ranking_scores = []
    with zipfile.ZipFile(zip_file, 'r') as fz:
        fz.extractall('temp_folder')
        for path in os.listdir('temp_folder'):
            long_path = os.path.join('temp_folder', path)
            name, ext = os.path.splitext(path)
            if "_summary_confidences_" in path and ext == ".json":
                with open(long_path, 'r') as file:
                    data = json.load(file)
                    ranking_score = data.get('ranking_score', None)
                    if ranking_score is not None:
                        ranking_scores.append({
                            'File Name': path,
                            'Ranking Score': ranking_score
                        })

    return pd.DataFrame(ranking_scores)

def extract_data_from_zip(zip_file):
    pdb_content = None
    ptmscore = None
    with zipfile.ZipFile(zip_file, 'r') as fz:
        fz.extractall('temp_folder')

    for path in os.listdir('temp_folder'):
        long_path = os.path.join('temp_folder', path)
        if long_path.endswith("_summary_confidences_0.json"):
            with open(long_path, 'r') as file:
                data = json.load(file)
                ptmscore = float(data['ptm'])

    for path in os.listdir('temp_folder'):
        long_path = os.path.join('temp_folder', path)
        if long_path.endswith("_model_0.cif"):
            file_parser = MMCIFParser(QUIET=True)
            structure = file_parser.get_structure("base", long_path)
            io = PDBIO()
            io.set_structure(structure)
            pdb_temp_file = "temp_folder/temp.pdb"
            io.save(pdb_temp_file)
            with open(pdb_temp_file, 'r') as file:
                pdb_content = file.read()

    for path in os.listdir('temp_folder'):
        long_path = os.path.join('temp_folder', path)
        if long_path.endswith("_full_data_0.json"):
            with open(long_path, 'r') as file:
                data = json.load(file)
                Pae_distance = pd.DataFrame(data['pae'])

    return pdb_content, Pae_distance, ptmscore

def visualize_pdb_3d(pdb_content):
    view_3d = py3Dmol.view(width=800, height=500)
    view_3d.addModel(pdb_content, 'pdb')
    view_3d.setStyle({'model': 0}, {"cartoon": {'color': '0x51adbe'}})
    view_3d.zoomTo()
    return view_3d

def GET_PAE_GRAPH(df_pae):
    fig_pae = go.Figure(data=go.Heatmap(z=df_pae, colorscale="Rainbow"))
    fig_pae.update_layout(width=500, height=500, yaxis=dict(scaleanchor="x", scaleratio=1), xaxis=dict(constrain='domain'))
    return fig_pae

def get_pdb_ca_from_content(pdb_content):
    Atom_serial_number = []
    Atom_name = []
    Residue_Name = []
    X_orthogonal_coordinates = []
    Y_orthogonal_coordinates = []
    Z_orthogonal_coordinates = []
    B_factor = []

    for linea in pdb_content.splitlines():
        if linea.startswith('ATOM') and linea[13:15].strip() == 'CA':
            Atom_serial_number.append(float(linea[6:11].strip()))
            Atom_name.append(linea[13:16].strip())
            Residue_Name.append(linea[17:20].strip())
            X_orthogonal_coordinates.append(float(linea[30:38].strip()))
            Y_orthogonal_coordinates.append(float(linea[38:46].strip()))
            Z_orthogonal_coordinates.append(float(linea[46:54].strip()))
            B_factor.append(float(linea[60:66].strip()))

    df_pdb = pd.DataFrame({
        'Atom Serial Number': Atom_serial_number,
        'Atom Name': Atom_name,
        'Residue Name': Residue_Name,
        'X orthogonal coordinate': X_orthogonal_coordinates,
        'Y orthogonal coordinate': Y_orthogonal_coordinates,
        'Z orthogonal coordinate': Z_orthogonal_coordinates,
        'pLDDT': B_factor
    })

    return df_pdb

def GET_PLDDTS(pdb_1):
    fig_plddt = go.Figure()
    fig_plddt.add_trace(go.Scatter(x=pdb_1.index, y=pdb_1["pLDDT"], line=dict(color='#40e0d0')))
    return fig_plddt

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

st.set_page_config(
    page_title="FoldXplore Tool",
    page_icon=":mag_right:",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# Título de la aplicación
st.title("FoldXplore")
st.write("-----")
st.markdown("""
    Bienvenido a FoldXplore! Aquí puedes explorar y visualizar estructuras de proteínas 
    obtenidas a partir de archivos de predicción de AlphaFold 3. Sube un archivo ZIP 
    y explora las predicciones de proteínas.
""")

#st.write("Bienvenidos a FoldXplore :)")

# Subir archivo ZIP
uploaded_file = st.file_uploader("Sube un archivo .zip de AlphaFold 3", type="zip")

if uploaded_file is not None:
    with open('temp_uploaded.zip', 'wb') as f:
        f.write(uploaded_file.getvalue())
    
    df_ranking_scores = extract_ranking_scores_from_zip('temp_uploaded.zip')
    
    if not df_ranking_scores.empty:
        st.write("Ranking Scores Extraídos:")
        st.dataframe(df_ranking_scores)
    else:
        st.write("No se encontraron archivos de resumen de confianza en el ZIP.")

    # Procesar archivo ZIP
    pdb_content, pae, ptmscore = extract_data_from_zip('temp_uploaded.zip')
    
    # Menú para seleccionar predicción
    prediction_option = st.selectbox("Selecciona la predicción a explorar", options=[0, 1, 2, 3, 4])

    if pdb_content:
        st.download_button(
            label="Descargar archivo PDB",
            data=pdb_content,
            file_name="model_relaxed.pdb",
            mime="chemical/x-pdb"
        )

    st.write('PTM Score: ', ptmscore)

    if pdb_content:
        st.write("Visualizando la estructura 3D del modelo PDB:")
        view_3d = visualize_pdb_3d(pdb_content)
        st.components.v1.html(view_3d._make_html(), height=500, width=800)

        st.write('PDB Dataframe')
        df_pdb = get_pdb_ca_from_content(pdb_content)
        st.write(df_pdb)
        st.write('PLDDT')
        fig_plddt = GET_PLDDTS(df_pdb)
        st.plotly_chart(fig_plddt)

        if pae is not None:
            st.write("PAE Heatmap:")
            fig_pae = GET_PAE_GRAPH(pae)
            st.plotly_chart(fig_pae)
    else:
        st.write("No se encontró un archivo _model_0.cif en el zip.")
