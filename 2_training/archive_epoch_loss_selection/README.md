# Archived: automatic epoch / loss selection from a sweep

Tools built (2026-11) to choose the training epoch and the loss of the MLPs from the `train_val_loss_curves.py` sweep (validation-loss minimum before overfitting, loss per
embedding and label by the MCC at that epoch), and a per-(embedding, label) epoch variant. **They are no longer part of the pipeline.** On the refreshed data they showed that

* the loss choice per label was within noise (winner margin 0.005 against a paired standard error of 0.009; leave-one-split-out it scored *below* plain BCE),
* the epoch chosen from a sweep on the held-out clusters is not independent of the test set, and different reasonable rules moved single embeddings by up to 0.02 MCC.

The pipeline now uses plain BCE and chooses the epoch inside every split by early stopping on a validation set that is cluster-disjoint from both train and test
(`../bakeoff.py`, see `../README.md`). The files are kept here, unchanged, because they document how those conclusions were reached; nothing imports them.
`test_select_epoch_*.py` run from this folder (`python test_select_epoch_and_loss.py`); the scripts read the output of `../train_val_loss_curves.py`.
