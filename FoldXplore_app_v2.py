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
    
    # Colorear por B-factor, asumiendo que los valores de B-factor están en el campo correspondiente
    # Azul (alto B-factor) -> Rojo (bajo B-factor)
    view_3d.setStyle({'model': 0}, {
        'cartoon': {
            'colorscheme': {
                'prop': 'b',  # Usar la B-factor como proxy para el color
                'gradient': 'roygb',  # Gradiente de rojo a azul
                'min': 0,  # Mínimo valor B-factor
                'max': 100  # Máximo valor B-factor
            }
        }
    })

    # Ajustar la visualización
    view_3d.zoomTo()
    return view_3d


def GET_PAE_GRAPH(df_pae):
    fig = px.imshow(df_pae, color_continuous_scale='Viridis', title="PAE Heatmap")
    
    fig.update_layout(
        width=500,
        height=500,
        xaxis_title="Scored Residue",  # Etiqueta del eje X
        yaxis_title="Aligned Residue",  # Etiqueta del eje Y
        coloraxis_colorbar=dict(
            title=dict(
                text="Expected Position Error (Ångströms)",  # Texto del título de la barra de colores
                side="bottom"  # Coloca el título debajo de la barra
            ),
            orientation="h",  # Orientación horizontal para la barra de colores
            x=0.5,  # Centrar la barra de colores debajo del gráfico
            y=-0.5,  # Posicionar la barra de colores más abajo
            len=1  # Extender la barra de colores para que coincida con el ancho del gráfico
        )
    )
    
    # Margen adicional para acomodar el título y la barra de colores
    fig.update_layout(margin=dict(l=50, r=50, t=50, b=120))
    
    return fig





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
    fig = px.line(pdb_1, x='Atom Serial Number', y="pLDDT", title="pLDDT Plot")
    fig.update_layout(
        xaxis_title="Cα positions ",  # Etiqueta del eje X
        yaxis_title="pLDDT"  # eje y 
    )
    return fig


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

    st.write('pTM Score: ', ptmscore)

    if pdb_content:
        st.write("Estructura 3D con colores basados en su estimacion de confianza pLDDT:")
        # Leyenda
        st.markdown("""
        <style>
            .legend {
                font-size: 16px;
                margin-bottom: 20px;
                background-color: #f9f9f9;
                padding: 10px;
                border-radius: 5px;
                border: 1px solid #ddd;
            }
        </style>
        <div class="legend">
            <strong>Leyenda de colores:</strong><br>
            <span style="color: blue;">&#9679;</span> Very high (pLDDT > 90)<br>
            <span style="color: lightblue;">&#9679;</span> Confident (90 > pLDDT > 70)<br>
            <span style="color: orange;">&#9679;</span> Low (70 > pLDDT > 50)<br>
            <span style="color: red;">&#9679;</span> Very low (pLDDT < 50)
        </div>
        """, unsafe_allow_html=True)
        view_3d = visualize_pdb_3d(pdb_content)
        # Renderizar el visor 3D en Streamlit
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