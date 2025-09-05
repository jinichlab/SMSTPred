from tabpfn import TabPFNClassifier
import os
import random
import subprocess
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.multioutput import MultiOutputClassifier
from sklearn.metrics import roc_auc_score, matthews_corrcoef
from xgboost import XGBClassifier
import torch
import torch.nn as nn
import torch.optim as optim
from tqdm import tqdm
from sklearn.svm import SVC
from sklearn.multioutput import MultiOutputClassifier
from sklearn.metrics import roc_auc_score, matthews_corrcoef, average_precision_score
from sklearn.metrics import precision_recall_curve, auc
from sklearn.model_selection import train_test_split
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.base import clone
from tqdm import tqdm
from joblib import Parallel, delayed
import psutil  # To check available memory
import faiss
import multiprocessing
import gc
from sklearn.base import clone

# Load the datasets
SAM_swiss_df = pd.read_csv('../files/SAM_swiss_substrate_types_balanced.csv')
# SAM_pubmed_df = pd.read_csv('../files/SAM_pubmed_enzymes_substrate_types.csv')
# UP_swiss_df = pd.read_csv('../files/UP_swiss_enzymes_substrate_types.csv')
# UP_pubmed_df = pd.read_csv('../files/UP_pubmed_enzymes_substrate_types.csv')

# 1) load entries to predict
predict_df = pd.read_excel("../files/uniprotkb_taxonomy_id_83332_AND_IPR0290_2025_06_23.xlsx")      # must have column "Entry"

# Print dataset shapes to verify
# print("SAM Swiss dataset shape:", SAM_swiss_df.shape)
# print("SAM PubMed dataset shape:", SAM_pubmed_df.shape)
# print("UP Swiss dataset shape:", UP_swiss_df.shape)
# print("UP PubMed dataset shape:", UP_pubmed_df.shape)

import h5py

def get_embeds(prott5_file, df_file):
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

# Define file paths for embeddings
SAM_swiss_prott5_file = '../files/SAM_swiss.h5'
# SAM_pubmed_prott5_file = '../files/SAM_swiss.h5'
# UP_swiss_prott5_file = '../files/per-protein.h5'
# UP_pubmed_prott5_file = '../files/per-protein.h5'
SAM_mtb_prott5_file = '../files/SAM_mtb.h5'

# Load embeddings into each dataframe
SAM_swiss_df = get_embeds(SAM_swiss_prott5_file, SAM_swiss_df)
# SAM_pubmed_df = get_embeds(SAM_pubmed_prott5_file, SAM_pubmed_df)
# UP_swiss_df = get_embeds(UP_swiss_prott5_file, UP_swiss_df)
# UP_pubmed_df = get_embeds(UP_pubmed_prott5_file, UP_pubmed_df)
predict_df = get_embeds(SAM_mtb_prott5_file, predict_df)

# Print updated dataset shapes to verify
# print("SAM Swiss dataset shape:", SAM_swiss_df.shape)
# print("SAM PubMed dataset shape:", SAM_pubmed_df.shape)
# print("UP Swiss dataset shape:", UP_swiss_df.shape)
# print("UP PubMed dataset shape:", UP_pubmed_df.shape)

datasets = {
    # "SAM_pubmed": SAM_pubmed_df,    
    "SAM_swiss": SAM_swiss_df,
    # "UP_pubmed": UP_pubmed_df,    
    # "UP_swiss": UP_swiss_df
}

# 🔹 Define classifiers with parallel processing
classifiers = {
    'XGBoost': XGBClassifier(eval_metric='logloss', random_state=42, n_jobs=-1),
    'Random Forest': RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1),
    'Logistic Regression': LogisticRegression(max_iter=1000, random_state=42, n_jobs=-1),
    'SVM_Linear': SVC(kernel='linear', probability=True, random_state=42),
    'SVM_RBF': SVC(kernel='rbf', probability=True, random_state=42),
    'SVM_Polynomial': SVC(kernel='poly', degree=3, probability=True, random_state=42),
    'SVM_Sigmoidal': SVC(kernel='sigmoid', probability=True, random_state=42),
    }

import pandas as pd
import numpy as np
from sklearn.multioutput import MultiOutputClassifier
from sklearn.base import clone

# ─── 1) Map dataset → classifier to use ──────────────────────────
predict_config = {
    # "UP_swiss":  "XGBoost",
    "SAM_swiss": "SVM_Polynomial",
    # "SAM_pubmed": "SVM_Linear",
    # add more dataset→classifier pairs as needed
}

# ─── 2) Load the list of entries & true labels ──────────────────
# predict.csv must have columns: Entry, Other, Protein, RNA, DNA (0/1 values)
predict_entries = set(predict_df["Entry"])

labels  = ['Other', 'Protein', 'RNA', 'DNA']
predictions = []
accuracies  = []

# ─── 3) Loop over each dataset/classifier pair ──────────────────
for ds_name, clf_name in predict_config.items():
    df = datasets[ds_name]
    # a) split train vs. held-out
    # train_df = df[~df['Entry'].isin(predict_entries)].reset_index(drop=True)
    # test_df  = df[df['Entry'].isin(predict_entries)].reset_index(drop=True)
    # if test_df.empty:
    #     print(f"⚠️  No held-out entries for {ds_name}, skipping.")
    #     continue
    test_df = predict_df
    train_df = df

    # # b) merge to get ground truth labels aligned
    # test_truth = pd.merge(
    #     test_df[['Entry']],
    #     predict_df[['Entry'] + labels],
    #     on='Entry', how='left'
    # ).reset_index(drop=True)

    # c) build numpy arrays
    X_train = np.vstack(train_df['Embeddings'].values).astype('float32')
    y_train = train_df[labels].values
    X_test  = np.vstack(test_df['Embeddings'].values).astype('float32')
    # y_true  = test_truth[labels].values
    # entries = test_truth['Entry'].tolist()
    entries = test_df['Entry'].tolist()

    # d) train the specified classifier
    base_clf = classifiers[clf_name]
    model    = MultiOutputClassifier(clone(base_clf))
    model.fit(X_train, y_train)

    # e) predict probabilities & binary calls
    prob_arrays = [est.predict_proba(X_test)[:,1] for est in model.estimators_]
    prob_matrix = np.vstack(prob_arrays).T       # shape = (n_test, n_labels)
    pred_matrix = (prob_matrix >= 0.5).astype(int)

    # f) record per-entry, per-label predictions
    for entry, probs, preds, in zip(entries, prob_matrix, pred_matrix):
        for i,label in enumerate(labels):
            predictions.append({
                "Entry":      entry,
                "Dataset":    ds_name,
                "Classifier": clf_name,
                "Label":      label,
                "Prob":       float(probs[i]),
                "Pred":       int(preds[i]),
                # "Truth":      int(truths[i])
            })

    # # g) compute accuracy per label for this dataset/classifier
    # for i,label in enumerate(labels):
    #     acc = float((pred_matrix[:,i] == y_true[:,i]).mean())
    #     accuracies.append({
    #         "Dataset":    ds_name,
    #         "Classifier": clf_name,
    #         "Label":      label,
    #         "Accuracy":   acc
    #     })

# ─── 4) Save results ─────────────────────────────────────────────
pred_df = pd.DataFrame(predictions)
pred_df.to_csv("mtb_SAM_predictions.csv", index=False)
print(f"✅ Saved mtb_SAM_predictions.csv ({len(pred_df)} rows)")

# acc_df = pd.DataFrame(accuracies)
# acc_df.to_csv("accuracy.csv", index=False)
# print(f"✅ Saved accuracy.csv ({len(acc_df)} rows)")
