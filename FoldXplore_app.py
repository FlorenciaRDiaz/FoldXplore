# ~ ~ ~ ~ ~ ~ ~ ~ LIBRERIAS ~ ~ ~ ~ ~ ~ ~
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

# ~ ~ ~ ~ ~ ~ ~ ~ FUNCIONES ~ ~ ~ ~ ~ ~ ~

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
            print(f"PDB file saved in: {pdb_file_path}")

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

def visualize_pdb_3d(pdb_content, spin=False):
    view_3d = py3Dmol.view(width=800, height=500)
    view_3d.addModel(pdb_content, 'pdb')

    # Colorear por B-factor
    view_3d.setStyle({'model': 0}, {
        'cartoon': {
            'colorscheme': {
                'prop': 'b',
                'gradient': 'roygb',
                'min': 0,
                'max': 100
            }
        }
    })

    # Controlar la rotación
    rotate = st.checkbox("Spin")
    if rotate:
        view_3d.spin(True)  # Activar la rotación
    else:
        view_3d.spin(False)  # Desactivar la rotación
    view_3d.zoomTo()
    return view_3d

def GET_PAE_GRAPH(df_pae):
    fig = px.imshow(df_pae, color_continuous_scale='Viridis', title="")

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
            x=0.5,  # Posiciona la barra de colores hacia el lado izquierdo
            y=-0.5,  # Alinea la barra de colores más cerca del gráfico
            len=0.5  # Reduce el ancho de la barra de colores
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
    fig = px.line(pdb_1, x='Atom Serial Number', y="pLDDT", title="")
    fig.update_layout(
        xaxis_title="Cα positions",
        yaxis_title="pLDDT"
    )
    return fig

# * * * * * * * * INTERFAZ DE USUARIO * * * * * * * *
st.set_page_config(
    page_title="FoldXplore Tool",
    page_icon=":mag_right:",
    layout="wide",
    initial_sidebar_state="collapsed"
)

#st.title("FoldXplore 3")
# Mostrar el logo
st.image("Logo.png",
         width = 300)
st.markdown("""
    Welcome to FoldXplore!,  Here you can explore and visualize protein structures
    obtained from AlphaFold 3 prediction files. Upload a ZIP
    file and explore the protein predictions.
""")
# Título grande fuera del file_uploader
st.write("### Upload a .zip file of AlphaFold 3: ")
# Objeto de carga de archivos sin título adicional
uploaded_file = st.file_uploader("", type="zip")


if uploaded_file is not None:
    with open('temp_uploaded.zip', 'wb') as f:
        f.write(uploaded_file.getvalue())

    df_ranking_scores = extract_ranking_scores_from_zip('temp_uploaded.zip')

    if not df_ranking_scores.empty:
        st.write("**Ranking Scores:**")
        st.dataframe(df_ranking_scores)
    else:
        st.write("No trusted summary files were found in the ZIP.")

    pdb_content, pae, ptmscore = extract_data_from_zip('temp_uploaded.zip')

    st.write("### Select the prediction to explore: ")
    prediction_option = st.selectbox("", options=[0, 1, 2, 3, 4])
    # Mensaje de descarga del PDB
    # st.markdown("### Download PDB File")
    # st.write("You can download the PDB file clicking on the button below:")
    # if pdb_content:
    #     st.download_button(
    #         label="Download PDB file",
    #         data=pdb_content,
    #         file_name="model_relaxed.pdb",
    #         mime="chemical/x-pdb"
    #     )
    # Estilo simple para el botón
    import base64

# Supón que pdb_content es tu contenido PDB

    if pdb_content:
        # Convertir pdb_content a base64
        b64 = base64.b64encode(pdb_content.encode()).decode()

        st.markdown("### Download PDB File")
        st.write("You can download the PDB file by clicking on the button below:")

        # Estilo simple para el botón
        st.markdown(
            """
            <style>
            .download-button {
                background-color: #A1C6E7;  /* Color celeste pastel */
                color: black;  /* Color del texto negro */
                padding: 10px 20px;
                text-align: center;
                border-radius: 5px;
                display: inline-block;
                margin: 10px 0;
                text-decoration: none;
                font-size: 16px;
            }
            </style>
            """, unsafe_allow_html=True
        )

        # Botón de descarga
        st.markdown(
            '<a href="data:chemical/x-pdb;base64,{data}" download="model_relaxed.pdb" class="download-button">Download PDB file</a>'.format(data=b64),
            unsafe_allow_html=True
        )
    else:
        st.error("No PDB content available to download.")
    
    st.markdown("### PDB Dataframe")
    df_pdb = get_pdb_ca_from_content(pdb_content)
    st.write(df_pdb)

    st.write('### pTM Score ')
    # Contenedor para mostrar el pTM Score dentro del cuadro
    st.markdown(f"""
    <div style="
        display: flex;
        align-items: center;
        background-color: #f8f9fa;
        border: 1px solid #ddd;
        border-radius: 8px;
        padding: 10px 20px;
        box-shadow: 2px 2px 8px rgba(0,0,0,0.1);
        width: 180px;
    ">
        <div style="
            background-color: #bae8e8;
            width: 6px;
            height: 100%;
            border-radius: 4px;
            margin-right: 10px;
        "></div>
        <div>
            <div style="font-size: 24px; font-weight: bold; color: #343a40;">{ptmscore}</div>
        </div>
    </div>
    """, unsafe_allow_html=True)

    if pdb_content:
        st.markdown("### 3D Structure")
        #st.write("3D structure with colors based on pLDDT confidence estimation:")
        spin = st.checkbox("Spin", value=True)  # Toggle for enabling spin
        # Add color-coded legend
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
            .legend-item {
                display: inline-block;
                margin-right: 15px;
            }
        </style>
        <div class="legend">
            <strong>Leyenda de colores:</strong><br>
            <div class="legend-item"><span style="color: blue;">&#9679;</span> Very high (pLDDT > 90)</div>
            <div class="legend-item"><span style="color: lightblue;">&#9679;</span> Confident (90 > pLDDT > 70)</div>
            <div class="legend-item"><span style="color: yellow;">&#9679;</span> Low (70 > pLDDT > 50)</div>
            <div class="legend-item"><span style="color: orange;">&#9679;</span> Very low (pLDDT < 50)</div>
        </div>
        """, unsafe_allow_html=True)

        view_3d = visualize_pdb_3d(pdb_content,spin=spin)
        # Renderizar el visor 3D en Streamlit
        st.components.v1.html(view_3d._make_html(), height=500, width=800)

        # st.write('PDB Dataframe')
        #df_pdb = get_pdb_ca_from_content(pdb_content)
        # st.write(df_pdb)
        st.write('### PLDDT Plot')
        fig_plddt = GET_PLDDTS(df_pdb)
        st.plotly_chart(fig_plddt)

        if pae is not None:
            st.write("### Predicted aligned error (PAE) ")
            fig_pae = GET_PAE_GRAPH(pae)
            st.plotly_chart(fig_pae, use_container_width=True)

    shutil.rmtree('temp_folder')
else:
    st.warning("No ZIP file has been uploaded.")

