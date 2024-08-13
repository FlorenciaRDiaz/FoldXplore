from AlphaFoldXplore_3 import alphafoldxplore as afx
from AlphaFoldXplore_3 import prediction_results 
from zipfile import ZipFile
import os
import re
import pandas as pd
import ipywidgets as widgets
from IPython.display import display
from google.colab import files

############### Funciones para cargar archivos [af3].zip --> [af3].afxt : 
import os
import shutil
import json
from zipfile import ZipFile
from Bio.PDB import MMCIFParser, PDBIO
from prediction_results import prediction_results
from google.colab import files  # Solo si usas Google Colab

# Función para manejar la carga del archivo en Google Colab
def upload_zip_colab():
    uploaded = files.upload()
    input_directory = 'input_af3'
    if not os.path.exists(input_directory):
        os.makedirs(input_directory)

    for filename in uploaded.keys():
        file_path = os.path.join(input_directory, filename)
        with open(file_path, 'wb') as f:
            f.write(uploaded[filename])

    # Devuelve la ruta del archivo cargado
    return file_path

# Función principal que integra la carga de archivo y el procesamiento
def load_af3_interactive_colab():
    # Solicitar al usuario que cargue el archivo .zip desde la computadora
    print("Por favor, carga tu archivo .zip desde tu computadora local:")
    location = upload_zip_colab()  # Usamos la función de carga

    # Solicitar al usuario que ingrese un nombre para la estructura
    name = input("Por favor, ingrese un nombre para la estructura: ")

    # Verificar si el archivo es un .zip válido
    if not location.endswith('.zip'):
        print("Error: El archivo proporcionado no es un archivo .zip.")
        return
    
    # Llamar a la función original load_af3 con los parámetros proporcionados
    Z = load_af3(name, location)
    
    return Z

# Función load_af3 (tu función original)
def load_af3(name, location):
    Z = {}
    
    if os.path.exists(name):
        print(f"Error: la carpeta con el nombre '{name}' ya existe en esta ubicación. Por favor, elige otro nombre o elimina/renombra la carpeta.")
        return
    
    folder = os.path.basename(location[:-4])
    extract_folder = f'AF3_files/{folder}'
    results_folder = f'AF3_files/{folder}_fx'
    os.makedirs(results_folder, exist_ok=True)
    
    with ZipFile(location, 'r') as fz:
        fz.extractall(extract_folder)  # Extraer el archivo .zip de AF3
    
    # Procesar los archivos extraídos
    for path in os.listdir(extract_folder):
        long_path = os.path.join(extract_folder, path)
        if long_path.endswith("_summary_confidences_0.json"):  # Obtener pTMscore
            with open(long_path, 'r') as file:
                data = json.load(file)
                ptmscore = float(data['ptm'])

        if long_path.endswith("_model_0.cif"):  # Convertir el CIF en PDB
            file_parser = MMCIFParser()
            structure = file_parser.get_structure("base", long_path)
            io = PDBIO()
            io.set_structure(structure)
            io.save(f'{results_folder}/{name}_relaxed.pdb')  # Guardar el archivo PDB

        if long_path.endswith("_full_data_0.json"):  # Obtener PAE
            with open(long_path, 'r') as file:
                data = json.load(file)
                distance = {"distance": data['pae']}
            with open(f'{results_folder}/{name}_pae.json', 'w', encoding='utf-8') as f:
                json.dump(distance, f, ensure_ascii=False)

    # Crear el reporte en la carpeta con el nombre proporcionado
    directory = f'{name}/{name}.zip'
    os.makedirs(f'{name}', exist_ok=True)
    with open(f'{name}/{name}_report.txt', 'w', encoding='utf-8') as file:
        file.write(name + '\n')
        file.write(directory + '\n')
        file.write("0" + '\n')
        file.write("no info" + '\n')
        file.write("pTMScore=" + str(ptmscore) + '\n')
        file.write("version=af3")
    
    # Mover y comprimir archivos
    os.system(f"mv '{results_folder}/{name}_pae.json' '{name}/{name}_pae.json'")
    os.system(f"mv '{results_folder}/{name}_relaxed.pdb' '{name}/{name}_relaxed.pdb'")
    os.system(f"zip -FSr '{name}.zip' '{name}'")
    shutil.rmtree(f'{name}')
    
    # Crear la entrada de predicción y guardarla en el diccionario
    prediction_entry = prediction_results(name, directory, "0", "no info", ptmscore)
    Z['p1'] = prediction_entry
    
    # Crear el archivo de lista
    os.makedirs(f'{name}', exist_ok=True)
    with open(f'{name}/{name}_list.txt', 'w', encoding='utf-8') as file:
        for result in list(Z.values()):
            file.write(f"{result.directory}\n")
    
    # Finalizar y mover el archivo .afxt
    os.system(f"mv '{name}.zip' '{name}/'")
    os.system(f"zip -FSr -D '{name}.zip' '{name}'")
    os.system(f"mv '{name}.zip' '{name}.afxt'")
    
    print(f"Guardado en tu computadora local. Nombre: \"{name}.afxt\"")
    
    return Z

############################################################################################## 
def get_pTMscore(results):
    """
    Función para extraer los valores de pTMscore y los nombres de los objetos prediction_results en el diccionario results.
    
    Parameters:
    results (dict): Diccionario que contiene objetos prediction_results.
    
    Returns:
    pd.DataFrame: DataFrame con los nombres de los archivos y sus correspondientes valores de pTMscore.
    """
    # Lista para almacenar los datos
    data = []
    
    # Recorrer cada entrada en el objeto results
    for key, result in results.items():
        # Acceder al valor de pTMscore
        ptmscore = result.ptmscore
        # Acceder al nombre del archivo
        name = result.name
        
        # Añadir el nombre del archivo, pTMscore y key al diccionario
        data.append({'key': key, 'Variant': name, 'pTMscore': ptmscore})
    
    # Crear un DataFrame con los datos
    df = pd.DataFrame(data)
    
    return df
