import streamlit as st
import requests
from rdkit import Chem
from rdkit.Chem import AllChem, MolToPDBBlock
import py3Dmol
import openai

# Configura tu clave de OpenAI
openai.api_key = st.secrets['api_key']

def generar_codigo_molecular(solicitud):
    prompt = """
    Escribe un script en Python que:
    -genera una estructura de la molécula que se pide
    sigue el siguiente ejemplo,
    # Importamos las librerías necesarias
    import requests
    from rdkit import Chem
    from rdkit.Chem import AllChem, MolToPDBBlock
    import py3Dmol

    cid = "10480"  # CID para el azul de metileno
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/record/SDF/?record_type=3d"
    response = requests.get(url)

    if response.status_code == 200:
        sdf_filename = "azul_de_metileno.sdf"
        with open(sdf_filename, "w") as file:
            file.write(response.text)
        print(f"Archivo SDF descargado y guardado como azul_de_metileno.sdf.")
    else:
        print(f"Error al descargar el archivo SDF: {response.status_code}")

    try:  
        suppl = Chem.SDMolSupplier(sdf_filename)
        mol = suppl[0]
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)

        pdb_block = MolToPDBBlock(mol)

        viewer = py3Dmol.view(width=400, height=400)
        viewer.addModel(pdb_block, 'pdb')
        viewer.setStyle({'stick':{}})
        viewer.zoomTo()
        viewer.show()
    except:
        print("Error al generar la estructura 3D")
    """ + solicitud + """
    Recuerda, solo debes escribir el código, sin comentarios ni explicaciones. No uses Markdown, ponlo como texto plano.
    """

    # Llamar a la API de OpenAI para generar el código
    response = openai.Completion.create(
        engine="gpt-4o-mini",  # O el modelo que prefieras
        prompt=prompt,
        max_tokens=800,
        temperature=0
    )
    
    return response.choices[0].text.strip()

# Función para mostrar la molécula generada
def mostrar_molecula(mol):
    pdb_block = MolToPDBBlock(mol)
    viewer = py3Dmol.view(width=400, height=400)
    viewer.addModel(pdb_block, 'pdb')
    viewer.setStyle({'stick': {}})
    viewer.zoomTo()
    viewer.show()

# Interfaz de usuario en Streamlit
st.title("Generador de Estructuras Moleculares")
nombre_molecula = st.text_input("Introduce el nombre de la molécula:", "aspirin")

if st.button("Generar estructura"):
    with st.spinner('Generando código y ejecutando...'):
        # Generar el código de la molécula
        solicitud_usuario = f"muestra la molecula de {nombre_molecula}"
        codigo_generado = generar_codigo_molecular(solicitud_usuario)
        
        # Ejecutar el código generado
        try:
            exec(codigo_generado)  # Ejecutar el código generado en un entorno seguro
            st.success(f"Estructura de {nombre_molecula} generada.")
        except Exception as e:
            st.error(f"Error al ejecutar el código: {e}")
