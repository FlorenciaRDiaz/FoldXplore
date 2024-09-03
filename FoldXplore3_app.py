import streamlit as st
import zipfile
import os
import json
import pandas as pd
from Bio.PDB import MMCIFParser, PDBIO, Select
import plotly.graph_objects as go
import py3Dmol
import io

# Funciones

def extract_data_from_zip(zip_file):
    pdb_contents = {}
    ptmscores = {}
    pae_data = {}

    with zipfile.ZipFile(zip_file, 'r') as fz:
        fz.extractall('temp_folder')

    for path in os.listdir('temp_folder'):
        long_path = os.path.join('temp_folder', path)
        name, ext = os.path.splitext(path)

        if ext == "_summary_confidences_0.json":
            with open(long_path, 'r') as file:
                data = json.load(file)
                ptmscores[name] = float(data['ptm'])

        if ext == "_model_0.cif":
            file_parser = MMCIFParser(QUIET=True)
            structure = file_parser.get_structure("base", long_path)
            io = PDBIO()
            io.set_structure(structure)
            pdb_content = io.get_structure().to_pdbstring()
            pdb_contents[name] = pdb_content

        if ext == "_full_data_0.json":
            with open(long_path, 'r') as file:
                data = json.load(file)
                pae_data[name] = pd.DataFrame(data['pae'])

    return pdb_contents, ptmscores, pae_data

def visualize_pdb_3d(pdb_content):
    view_3d = py3Dmol.view(width=800, height=500)
    view_3d.addModel(pdb_content, 'pdb')
    view_3d.setStyle({'model': 0}, {"cartoon": {'color': '0x51adbe'}})
    view_3d.zoomTo()
    return view_3d

def GET_PAE_GRAPH(df_pae):
    fig_pae = go.Figure(data=go.Heatmap(z=df_pae, colorscale="Rainbow"))
    fig_pae.update_layout(
        width=500,
        height=500,
        yaxis=dict(scaleanchor="x", scaleratio=1),
        xaxis=dict(constrain='domain')
    )
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
        'B factor': B_factor
    })

    return df_pdb

def GET_PLDDTS(df_pdb):
    fig_plddt = go.Figure()
    fig_plddt.add_trace(go.Scatter(x=df_pdb.index, y=df_pdb["B factor"], line=dict(color='#40e0d0')))
    return fig_plddt

# Configuración de la aplicación
st.set_page_config(
    page_title="FoldXplore",
    page_icon=":mag_right:",
    layout="wide",
    initial_sidebar_state="collapsed"
)

st.write("Bienvenidos a FoldXplore :)")

uploaded_file = st.file_uploader("Sube un archivo .zip de AlphaFold 3", type="zip")

if uploaded_file is not None:
    with open('temp_uploaded.zip', 'wb') as f:
        f.write(uploaded_file.getvalue())
    
    pdb_contents, ptmscores, pae_data = extract_data_from_zip('temp_uploaded.zip')
    
    if pdb_contents:
        selected_seed = st.selectbox("Selecciona la semilla de predicción", options=list(pdb_contents.keys()))

        if selected_seed:
            st.write(f"Seleccionaste la semilla: {selected_seed}")

            st.write('PTM Score: ', ptmscores.get(selected_seed, "No disponible"))

            st.write("Visualizando la estructura 3D del modelo PDB:")
            view_3d = visualize_pdb_3d(pdb_contents[selected_seed])
            st.components.v1.html(view_3d._make_html(), height=500, width=800)

            st.write('PDB Dataframe')
            df_pdb = get_pdb_ca_from_content(pdb_contents[selected_seed])
            st.write(df_pdb)

            st.write('PLDDT')
            fig_plddt = GET_PLDDTS(df_pdb)
            st.plotly_chart(fig_plddt)

            if selected_seed in pae_data:
                st.write('PAE Heatmap:')
                fig_pae = GET_PAE_GRAPH(pae_data[selected_seed])
                st.plotly_chart(fig_pae)
            else:
                st.write("No se encontró datos PAE para la semilla seleccionada.")

            st.download_button(
                label="Descargar archivo PDB",
                data=pdb_contents[selected_seed],
                file_name=f"{selected_seed}_model_relaxed.pdb",
                mime="chemical/x-pdb"
            )
    else:
        st.write("No se encontraron archivos PDB en el ZIP.")
