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
from Bio.PDB import MMCIFParser, PDBIO
from Bio.PDB import MMCIFParser, PDBIO, Select
from io import BytesIO

import py3Dmol


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 

# * * * * * * * * FUNCIONES * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  

# FUNCION- FLOR -----------------------------------------------------------------------------------------------
def load_af3_pdb(zip_file_path):
    Z = {}
    name = os.path.basename(zip_file_path).replace('.zip', '')
    folder = f'AF3_files/{name}'
    results_folder = f'AF3_files/{name}'
    os.makedirs(results_folder, exist_ok=True)

    with zipfile.ZipFile(zip_file_path, 'r') as fz:
        fz.extractall(folder)

    pdb_file_path = None

    for path in os.listdir(folder):
        long_path = os.path.join(folder, path)

        if long_path.endswith("_summary_confidences_0.json"):
            with open(long_path, 'r') as file:
                data = json.load(file)
                ptmscore = float(data['ptm'])

        if long_path.endswith("_model_0.cif"):
            file_parser = MMCIFParser(QUIET=True)
            structure = file_parser.get_structure("base", long_path)
            io = PDBIO()
            io.set_structure(structure)
            io.save(f'{results_folder}/{name}_relaxed.pdb', select=Select())
            pdb_file_path = f'{results_folder}/{name}_relaxed.pdb'
            print(f"Archivo PDB guardado en: {pdb_file_path}")

        if long_path.endswith("_full_data_0.json"):
            with open(long_path, 'r') as file:
                data = json.load(file)
                distance = {"distance": data['pae']}
            with open(f'{results_folder}/{name}_pae.json', 'w', encoding='utf-8') as f:
                json.dump(distance, f, ensure_ascii=False)

    return pdb_file_path
# def load_af3_pdb(zip_file_path):  # Cambié el argumento para aceptar la ruta del archivo
#     Z = {}

#     # Obtener el nombre del archivo sin extensión para usarlo como nombre de la carpeta
#     name = os.path.basename(zip_file_path).replace('.zip', '')

#     folder = f'AF3_files/{name}'
#     results_folder = f'AF3_files/{name}'
#     os.makedirs(results_folder, exist_ok=True)

#     # Extraer el archivo ZIP
#     with zipfile.ZipFile(zip_file_path, 'r') as fz:
#         fz.extractall(folder)

#     pdb_file_path = None  # Inicializar la ruta del archivo PDB

#     # Procesar los archivos dentro de la carpeta extraída
#     for path in os.listdir(folder):
#         long_path = os.path.join(folder, path)

#         if long_path.endswith("_summary_confidences_0.json"):
#             with open(long_path, 'r') as file:
#                 data = json.load(file)
#                 ptmscore = float(data['ptm'])

#         if long_path.endswith("_model_0.cif"):
#             file_parser = MMCIFParser(QUIET=True)
#             structure = file_parser.get_structure("base", long_path)
#             io = PDBIO()
#             io.set_structure(structure)
#             io.save(f'{results_folder}/{name}_relaxed.pdb', select=Select())  # Asegurar que select esté configurado
#             pdb_file_path = f'{results_folder}/{name}_relaxed.pdb'
#             print(f"Archivo PDB guardado en: {pdb_file_path}")

#         if long_path.endswith("_full_data_0.json"):
#             with open(long_path, 'r') as file:
#                 data = json.load(file)
#                 distance = {"distance": data['pae']}
#             with open(f'{results_folder}/{name}_pae.json', 'w', encoding='utf-8') as f:
#                 json.dump(distance, f, ensure_ascii=False)

#     # Retornar la ruta del archivo PDB
#     return pdb_file_path

#--------------------------------------------------------------------------------------------------------------

def extract_data_from_zip(zip_file):


    pdb_content = None  # Para almacenar el contenido PDB
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

            # Guardar el archivo PDB en un archivo temporal en disco
            pdb_temp_file = "temp_folder/temp.pdb"
            io.save(pdb_temp_file)  # Guardar en archivo temporal

            # Leer el contenido del archivo temporal en memoria
            with open(pdb_temp_file, 'r') as file:
                pdb_content = file.read()

    
    
    for path in os.listdir('temp_folder'):
        long_path = os.path.join('temp_folder', path)        
        
        if long_path.endswith("_full_data_0.json"):
            with open(long_path, 'r') as file:
                data = json.load(file)
                Pae_distance = pd.DataFrame(data['pae'])
                
            
    return pdb_content, Pae_distance, ptmscore


# Función para visualizar el archivo PDB en 3D
def visualize_pdb_3d(pdb_content):
    view_3d = py3Dmol.view(width=800, height=500)
    view_3d.addModel(pdb_content, 'pdb')
    view_3d.setStyle({'model': 0}, {"cartoon": {'color': '0x51adbe'}})  # Modelo 0 con color azul
    view_3d.zoomTo()
    #view_3d.show()
    return view_3d

def GET_PAE_GRAPH(df_pae):
    fig_pae = go.Figure(data=go.Heatmap(
                    z=df_pae, 
                    colorscale="Rainbow"))
    # Obtener la cantidad de filas y columnas del dataframe
    # num_rows, num_cols = df_pae.shape
    
    # Ajustar el layout para hacer el heatmap cuadrado
    fig_pae.update_layout(
        width=500,  # Establecer el ancho deseado
        height=500,  # Establecer la altura deseada
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
            Atom_name.append(linea[13:16].strip())  # Strip elimina los espacios en blanco en el string
            Residue_Name.append(linea[17:20].strip())
            X_orthogonal_coordinates.append(float(linea[30:38].strip()))
            Y_orthogonal_coordinates.append(float(linea[38:46].strip()))
            Z_orthogonal_coordinates.append(float(linea[46:54].strip()))
            B_factor.append(float(linea[60:66].strip()))

    # DATAFRAME CON LOS VALORES DEL PDB- CA
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
    fig_plddt.add_trace(go.Scatter( x=pdb_1.index,
                                    y=pdb_1["pLDDT"],
                                    line=dict(color='#40e0d0')))

    return fig_plddt

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 


# Set Up de la applicación * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
st.set_page_config(
    page_title="FoldXplore",
    page_icon="	:mag_right:",
    layout="wide",
    initial_sidebar_state="collapsed",
    menu_items={
        # #pestañas del menu que se pueden modificar, pero solo acepta estas 3
        # 'Get Help': "COMPLETAR MAIL",
        # 'Report a bug': "COMPLETAR - MAIL",
        # 'About': "COMPLETAR"
    }
)
#objeto predeterminaod logo de pestaña, 
#   FUNCION GET_PLDDTS
#Parametros
    #pdb1 = dataframe con el primer pdb
    #ref_1 = string con referencia al pdf anterior
    #pdb2 = dataframe con el segundo pdb
    #ref_2 = string con referencia al pdf anterior
#Return
    #fig_plddt




# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

st.write("Bienvenidos a FoldXplore :)")


# objeto para subir el archivo .zp
uploaded_file = st.file_uploader("Sube un archivo .zip de AlphaFold 3", type="zip")

if uploaded_file is not None:
     # Guardar el archivo subido temporalmente
    with open('temp_uploaded.zip', 'wb') as f:
        f.write(uploaded_file.getvalue())
    
   #procesar el archivo ZIP y extraer los datos necesarios
    pdb_content, pae, ptmscore = extract_data_from_zip(uploaded_file)
    
    
    # Botón de descarga para el archivo PDB
    if pdb_content:
        st.download_button(
            label="Descargar archivo PDB",
            data=pdb_content,
            file_name="model_relaxed.pdb",
            mime="chemical/x-pdb"
        )

    st.write('PTM Score: ', ptmscore)
    
    st.write("PAE's Dataframe")
    st.dataframe(pae)

    fig_pae = GET_PAE_GRAPH(pae)
    st.write('PAE´s Heatmap:')
    st.plotly_chart(fig_pae)

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
    else:
        st.write("No se encontró un archivo _model_0.cif en el zip.")
