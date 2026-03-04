import pandas as pd
import h5py
import numpy as np

def get_embeds (prott5_file, df_file):
    entries_list = []
    embeddings_list = []
    data_map = {}

    with h5py.File(prott5_file, "r") as file:
        for sequence_id, embedding in file.items():
            entries_list.append(sequence_id)
            embeddings = np.array(embedding)
            embeddings_list.append(embeddings)
            data_map[sequence_id] = embeddings

    prott5_df = df_file
    prott5_df['Embeddings'] = prott5_df['Entry'].map(data_map)
    prott5_df['Embeddings'] = prott5_df['Embeddings'].apply(lambda x: x.tolist())
    #Remove proteins without embeddings from the dataframe
    prott5_df = prott5_df.dropna(subset=['Embeddings'])

    return(prott5_df)

df = pd.read_csv('SAM_swiss_small_molecules_smiles.csv')
prott5_file = 'example.h5'
df = get_embeds(prott5_file, df)

file_out = 'SAM_Swiss_SM_smiles_Embeddings.csv'

df.to_csv(file_out, index = False)
