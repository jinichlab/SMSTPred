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

# Load the datasets
# SAM_swiss_df = pd.read_csv('../files/SAM_swiss_enzymes_substrate_types.csv')
# SAM_pubmed_df = pd.read_csv('../files/SAM_pubmed_enzymes_substrate_types.csv')
UP_swiss_df = pd.read_csv('../files/UP_swiss_enzymes_substrate_types.csv')
# UP_pubmed_df = pd.read_csv('../files/UP_pubmed_enzymes_substrate_types.csv')

# Print dataset shapes to verify
# print("SAM Swiss dataset shape:", SAM_swiss_df.shape)
# print("SAM PubMed dataset shape:", SAM_pubmed_df.shape)
print("UP Swiss dataset shape:", UP_swiss_df.shape)
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
# SAM_swiss_prott5_file = '../files/SAM_swiss.h5'
# SAM_pubmed_prott5_file = '../files/SAM_swiss.h5'
UP_swiss_prott5_file = '../files/per-protein.h5'
# UP_pubmed_prott5_file = '../files/per-protein.h5'

# Load embeddings into each dataframe
# SAM_swiss_df = get_embeds(SAM_swiss_prott5_file, SAM_swiss_df)
# SAM_pubmed_df = get_embeds(SAM_pubmed_prott5_file, SAM_pubmed_df)
UP_swiss_df = get_embeds(UP_swiss_prott5_file, UP_swiss_df)
# UP_pubmed_df = get_embeds(UP_pubmed_prott5_file, UP_pubmed_df)

# Print updated dataset shapes to verify
# print("SAM Swiss dataset shape:", SAM_swiss_df.shape)
# print("SAM PubMed dataset shape:", SAM_pubmed_df.shape)
print("UP Swiss dataset shape:", UP_swiss_df.shape)
# print("UP PubMed dataset shape:", UP_pubmed_df.shape)

datasets = {
    # "SAM_pubmed": SAM_pubmed_df,    
    # "SAM_swiss": SAM_swiss_df,
    # "UP_pubmed": UP_pubmed_df,    
    "UP_swiss": UP_swiss_df
}

# 🔹 Ensure PyTorch uses 20 CPU threads
torch.set_num_threads(20)
device = torch.device("cpu")

# 🔹 Function to estimate memory usage
def estimate_memory_usage(n_samples):
    """Estimate memory required for a full cosine similarity matrix."""
    memory_needed_gb = (n_samples ** 2 * 4) / (1024 ** 3)  # Convert bytes to GB (float32)
    return memory_needed_gb

# 🔹 Compute cosine similarity in batches & store in disk-backed array
def batched_cosine_similarity(embeddings, batch_size=5000, save_path="cosine_sim_matrix.dat"):
    """Compute cosine similarity in batches & store it in a disk-backed file (memory-efficient)."""
    n = embeddings.shape[0]
    sim_matrix = np.memmap(save_path, dtype=np.float32, mode='w+', shape=(n, n))  # ✅ Disk-backed storage

    for i in tqdm(range(0, n, batch_size), desc="🔹 Computing Cosine Similarity (Batched)", leave=True):
        end_i = min(i + batch_size, n)
        for j in range(0, n, batch_size):
            end_j = min(j + batch_size, n)
            sim_matrix[i:end_i, j:end_j] = cosine_similarity(embeddings[i:end_i], embeddings[j:end_j])

    return sim_matrix

# 🔹 Compute cosine similarity on demand (Greedy, Sorted by Norm)
def on_demand_cosine_similarity(embeddings, threshold):
    """Greedy clustering with on-demand similarity calculation, sorting by L2 norm."""
    n = len(embeddings)
    embedding_norms = np.linalg.norm(embeddings, axis=1)
    sorted_indices = np.argsort(-embedding_norms)  # Sort by descending norm

    clusters = {}
    cluster_labels = -1 * np.ones(n, dtype=int)
    cluster_id = 0

    # Iterate over sorted embeddings
    for i in tqdm(sorted_indices, desc="On-Demand Greedy Clustering", leave=True):
        if cluster_labels[i] != -1:
            continue  # Already assigned

        # Start a new cluster with this embedding as the representative
        clusters[cluster_id] = [i]
        cluster_labels[i] = cluster_id  

        # Compare only to existing cluster representatives
        for rep in clusters.keys():
            rep_index = clusters[rep][0]  # First element of each cluster is the representative
            sim = cosine_similarity([embeddings[i]], [embeddings[rep_index]])[0, 0]

            if sim >= threshold:
                cluster_labels[i] = rep
                clusters[rep].append(i)
                break  # Stop checking once assigned

        if cluster_labels[i] == cluster_id:  # If still unassigned, it stays as a new cluster
            cluster_id += 1  

    return clusters, cluster_labels


# 🔹 Greedy clustering from precomputed similarity matrix (Sorted by Norm)
def greedy_clustering_from_matrix(sim_matrix, threshold):
    """Greedy clustering using a **precomputed** cosine similarity matrix, sorting by L2 norm."""
    n = len(sim_matrix)
    embedding_norms = np.linalg.norm(sim_matrix, axis=1)  # Estimate norms from similarities
    sorted_indices = np.argsort(-embedding_norms)  # Sort in descending order

    clusters = {}
    cluster_labels = -1 * np.ones(n, dtype=int)
    cluster_id = 0

    # Iterate over sorted embeddings
    for i in tqdm(sorted_indices, desc="Greedy Clustering (Precomputed)", leave=True):
        if cluster_labels[i] != -1:
            continue  # Already assigned

        # Start a new cluster with this embedding as the representative
        clusters[cluster_id] = [i]
        cluster_labels[i] = cluster_id  

        # Compare only to existing cluster representatives
        for rep in clusters.keys():
            rep_index = clusters[rep][0]  # First element of each cluster is the representative
            sim = sim_matrix[i, rep_index]  # Retrieve similarity from precomputed matrix

            if sim >= threshold:
                cluster_labels[i] = rep
                clusters[rep].append(i)
                break  # Stop checking once assigned

        if cluster_labels[i] == cluster_id:  # If still unassigned, it stays as a new cluster
            cluster_id += 1  

    return clusters, cluster_labels

# 🔹 Function to determine clustering strategy (stored matrix vs. on-demand)
def dynamic_clustering(embeddings, threshold):
    """Automatically determines whether to store the full similarity matrix or compute on demand."""
    n_samples = len(embeddings)
    available_memory_gb = psutil.virtual_memory().available / (1024 ** 3)  # ✅ Convert to GB
    required_memory_gb = estimate_memory_usage(n_samples)

    print(f"\nDataset Size: {n_samples} samples")
    print(f"Memory Needed for Full Matrix: {required_memory_gb:.2f} GB")
    print(f"Available System Memory: {available_memory_gb:.2f} GB")

    if required_memory_gb < available_memory_gb * 0.5:  # ✅ Use full similarity matrix if memory allows
        print("Using Precomputed Cosine Similarity Matrix (Fast Clustering)")
        sim_matrix = batched_cosine_similarity(embeddings, 5000)
        return greedy_clustering_from_matrix(sim_matrix, threshold)
    
    else:  # ✅ Use on-demand similarity if memory is limited
        print("Memory Constraint! Switching to On-Demand Similarity Calculation")
        return on_demand_cosine_similarity(embeddings, threshold)

# 🔹 Function to split clusters into train and test sets
def split_clusters(clusters, train_ratio=0.8):
    """
    Splits the clusters into training and testing sets based on a ratio.
    """
    cluster_ids = list(clusters.keys())
    random.shuffle(cluster_ids)

    train_size = int(train_ratio * len(cluster_ids))
    train_clusters = {k: clusters[k] for k in cluster_ids[:train_size]}
    test_clusters = {k: clusters[k] for k in cluster_ids[train_size:]}

    train_indices = [idx for cluster in train_clusters.values() for idx in cluster]
    test_indices = [idx for cluster in test_clusters.values() for idx in cluster]

    return train_indices, test_indices

# 🔹 Define classifiers with parallel processing
classifiers = {
    'XGBoost': XGBClassifier(eval_metric='logloss', random_state=42, n_jobs=-1),
    'Random Forest': RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1),
    'Logistic Regression': LogisticRegression(max_iter=1000, random_state=42, n_jobs=-1),
    # 'SVM_Linear': SVC(kernel='linear', probability=True, random_state=42),
    # 'SVM_RBF': SVC(kernel='rbf', probability=True, random_state=42),
    # 'SVM_Polynomial': SVC(kernel='poly', degree=3, probability=True, random_state=42),
    # 'SVM_Sigmoidal': SVC(kernel='sigmoid', probability=True, random_state=42),
    }

# 🔹 Neural Network Model
class MultilabelNN(nn.Module):
    def __init__(self, input_size, output_size):
        super(MultilabelNN, self).__init__()
        self.fc1 = nn.Linear(input_size, 512)
        self.relu = nn.ReLU()
        self.fc2 = nn.Linear(512, 256)
        self.fc3 = nn.Linear(256, output_size)
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        x = self.fc1(x)
        x = self.relu(x)
        x = self.fc2(x)
        x = self.relu(x)
        x = self.fc3(x)
        return self.sigmoid(x)

def train_svm_parallel(X_train, y_train, classifier_name):
    """
    Train multiple SVM classifiers in parallel using the predefined classifiers dictionary.
    """
    # Retrieve the predefined SVM model
    if classifier_name not in classifiers:
        raise ValueError(f"Unknown SVM classifier: {classifier_name}")

    base_svm = classifiers[classifier_name]  # Predefined SVM model from `classifiers`
    
    def train_single_svm(label_idx):
        if len(np.unique(y_train[:, label_idx])) == 1:  # Skip training if only one class
            return None  
        model = clone(base_svm)  # Clone the predefined SVM model
        model.fit(X_train, y_train[:, label_idx])
        return model

    trained_svms = Parallel(n_jobs=-1)(delayed(train_single_svm)(i) for i in range(y_train.shape[1]))
    return [svm for svm in trained_svms if svm is not None]  # Remove None values

def evaluate_classifier(df, df_name, classifier, classifier_name, threshold, clusters, iterations=5, epochs=1000):
    """
    Evaluates a classifier by training and testing it across multiple iterations.
    Now uses precomputed clusters to avoid recomputation.
    """
    labels = ['Other', 'Protein', 'RNA', 'DNA']
    metrics = []

    embeddings = np.vstack(df['Embeddings'].values).astype(np.float32)
    for i in tqdm(range(iterations), desc=f"Running {classifier_name} on {df_name} at {threshold:.2f}", leave=True):

        # Try splitting until all labels have both 0 and 1
        max_retries = 10
        for attempt in range(max_retries):
            train_indices, test_indices = split_clusters(clusters, train_ratio=0.8)
            df_train_tmp = df.iloc[train_indices]

            valid_split = True
            for label in labels:
                unique_vals = np.unique(df_train_tmp[label])
                if len(unique_vals) < 2:
                    valid_split = False
                    break  # No need to check further if any label is invalid

            if valid_split:
                break  # Good split found, exit retry loop
        else:
            print(f"Skipping iteration {i+1} for {df_name} at threshold {threshold:.2f}: could not get a valid train split after {max_retries} attempts.")
            continue  # After retries, if still bad, skip iteration

        # ✅ Check if test set is empty
        if len(test_indices) == 0:
            print(f"Skipping {classifier_name} on {df_name} at {threshold:.2f} due to empty test set.")
            continue  # ✅ Skip this iteration

        df_train = df.iloc[train_indices].reset_index(drop=True)
        df_test = df.iloc[test_indices].reset_index(drop=True)

        X_train = np.vstack(df_train['Embeddings'].values)
        y_train = df_train[labels].values
        X_test = np.vstack(df_test['Embeddings'].values)
        y_test = df_test[labels].values

        if classifier_name == 'Neural Network':
            
            # ------------------------------------------------------------------
            # carve training fold into (train_inner ,  val)  – 10 % is typical
            # ------------------------------------------------------------------
            X_train_nn, X_val_nn, y_train_nn, y_val_nn = train_test_split(
                    X_train, y_train, test_size=0.10, random_state=42, stratify=None)

            # tensors
            X_train_tensor = torch.tensor(X_train_nn, dtype=torch.float32)
            y_train_tensor = torch.tensor(y_train_nn, dtype=torch.float32)
            X_val_tensor   = torch.tensor(X_val_nn,   dtype=torch.float32)
            y_val_tensor   = torch.tensor(y_val_nn,   dtype=torch.float32)

            X_test_tensor  = torch.tensor(X_test,     dtype=torch.float32)   # untouched
            y_test_tensor  = torch.tensor(y_test,     dtype=torch.float32)


            # ✅ Initialize the Neural Network
            model = MultilabelNN(input_size=1024, output_size=len(labels))
            criterion = nn.BCELoss()
            optimizer = optim.Adam(model.parameters(), lr=0.001)

            best_epoch = -1
            best_model_path = "best_model_path.pth"
            epochs = 100
            best_auprc  = -np.inf
            best_val_loss = float('inf')
            epochs_no_improve = 0
            patience = 5          # stop after 5 non-improving epochs

            for epoch in range(epochs):
                # ---- Train step ----
                model.train()
                optimizer.zero_grad()
                outputs_train = model(X_train_tensor)
                loss_train = criterion(outputs_train, y_train_tensor)
                loss_train.backward()
                optimizer.step()

                # ---- Validation step ----
                model.eval()
                with torch.no_grad():
                    outputs_val = model(X_val_tensor)
                    val_loss = criterion(outputs_val, y_val_tensor).item()
                    y_val_probs = outputs_val.numpy()
                    prec, rec, _ = precision_recall_curve(y_val_tensor.numpy().flatten(),
                                                          y_val_probs.flatten())
                    auprc = auc(rec, prec)

                # ---- Checkpoint / early-stop ----
                if val_loss < best_val_loss:
                    best_val_loss = val_loss
                    best_auprc = auprc
                    best_epoch = epoch
                    torch.save(model.state_dict(), best_model_path)
                    epochs_no_improve = 0
                    print(f"✅  epoch {epoch+1}: new best val-loss {val_loss:.4f}")
                else:
                    epochs_no_improve += 1
                    if epochs_no_improve >= patience:
                        print(f"⏹️  early stopping at epoch {epoch+1} (no val improvement)")
                        break

                print(f"Epoch {epoch+1:3d} | train-loss {loss_train.item():.4f} | "
                      f"val-loss {val_loss:.4f} | val-AUPRC {auprc:.4f}")

            best_model = MultilabelNN(input_size=1024, output_size=len(labels))
            best_model.load_state_dict(torch.load(best_model_path))
            best_model.eval()

            with torch.no_grad():
                y_pred_probs = best_model(X_test_tensor).numpy()
                y_pred = (y_pred_probs >= 0.5).astype(int)

        elif classifier_name.startswith("SVM"):
            trained_svms = train_svm_parallel(X_train, y_train, classifier_name)  # ✅ Use predefined classifiers
            if len(trained_svms) == 0:
                print(f"Warning: No valid SVM models trained for {classifier_name} on {df_name} at threshold {threshold}")
                y_pred_probs = np.zeros((X_test.shape[0], len(labels)))
            else:
                y_pred_probs = np.array([svm.predict_proba(X_test)[:, 1] for svm in trained_svms]).T

            if y_pred_probs.shape[1] != len(labels):
                print(f"Shape Mismatch: Expected {len(labels)} labels, but got {y_pred_probs.shape[1]}")
                missing_columns = len(labels) - y_pred_probs.shape[1]
                y_pred_probs = np.hstack([y_pred_probs, np.zeros((X_test.shape[0], missing_columns))])
            y_pred = (y_pred_probs >= 0.5).astype(int)

        else:
            multi_output_model = MultiOutputClassifier(classifier)
            multi_output_model.fit(X_train, y_train)
            y_pred = multi_output_model.predict(X_test)
            y_prob_list = [est.predict_proba(X_test) for est in multi_output_model.estimators_]
            y_pred_probs = np.hstack([prob[:, 1].reshape(-1, 1) for prob in y_prob_list])

        for idx, label in enumerate(labels):
            y_test_binary = y_test[:, idx]
            y_pred_binary = y_pred[:, idx]
            y_prob_binary = y_pred_probs[:, idx]

            # Compute AUROC
            auroc = roc_auc_score(y_test_binary, y_prob_binary) if len(np.unique(y_test_binary)) > 1 else np.nan

            # Compute MCC
            mcc = matthews_corrcoef(y_test_binary, y_pred_binary) if len(np.unique(y_test_binary)) > 1 else np.nan

            # Compute AP (Average Precision Score)
            ap = average_precision_score(y_test[:, idx], y_pred_probs[:, idx])

            # Compute AUPRC (Area Under Precision-Recall Curve) explicitly
            precision, recall, _ = precision_recall_curve(y_test[:, idx], y_pred_probs[:, idx])
            auprc = auc(recall, precision)

            # Append all metrics
            metrics.append({
                'Dataset': df_name,
                'Classifier': classifier_name,
                'Threshold': threshold,
                'Label': label,
                'Iteration': i + 1,
                'AUROC': auroc,
                'MCC': mcc,
                'AP': ap,
                'AUPRC': auprc
            })

    return metrics

thresholds = np.arange(0.1, 1.05, 0.1)

file_path = "classifier_results_all_datasets.csv"

# 🔹 Ensure old file is removed before writing new results
if os.path.exists(file_path):
    os.remove(file_path)

results_df = pd.DataFrame()  # Initialize an empty DataFrame

for df_name, df in datasets.items():
    embeddings = np.vstack(df['Embeddings'].values).astype(np.float32)

    for threshold in tqdm(thresholds, desc=f"Processing {df_name}"):
        clusters, cluster_labels = dynamic_clustering(embeddings, threshold)

        threshold_results = []  # Collect results for this threshold only

        for name, clf in classifiers.items():
            new_results = evaluate_classifier(df, df_name, clf, name, threshold, clusters)
            threshold_results.extend(new_results)

        # Run Neural Network separately
        new_results_nn = evaluate_classifier(df, df_name, None, 'Neural Network', threshold, clusters)
        threshold_results.extend(new_results_nn)

        # Convert and write only current threshold's results
        threshold_df = pd.DataFrame(threshold_results)
        threshold_df.to_csv(file_path, mode='a', header=not os.path.exists(file_path), index=False)

        # Also keep cumulative results in memory if needed
        results_df = pd.concat([results_df, threshold_df], ignore_index=True)

        # ✅ Cleanup to save memory
        import gc
        for var in ['sim_matrix', 'cluster_labels', 'clusters']:
            if var in locals():
                del locals()[var]
        gc.collect()