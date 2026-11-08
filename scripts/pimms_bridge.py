"""Bridge to PIMMS (pimms-learn 0.5.2; Webel et al. 2024, doi:10.1038/s41467-024-48711-5), run in its own virtualenv so
its torch/fastai stack never touches the benchmark environment.

    <pimms-venv>/bin/python scripts/pimms_bridge.py in.npz out.npz {cf|vae} seed

in.npz: X (cells x proteins, log2, NaN = missing). out.npz: Y (full matrix). Settings follow the PIMMS README/tutorial
defaults: CF n_factors=30, epochs_max=20, early stopping (patience 1) on a 5% validation split of observed entries;
VAE hidden_layers=[64], latent_dim=10, epochs_max=300, patience 50 on the same kind of split. CPU only.
"""

import os
import sys
import tempfile
import warnings

warnings.filterwarnings("ignore")
os.environ.setdefault("OMP_NUM_THREADS", "1")


def main():
    import numpy as np
    import pandas as pd
    import torch

    torch.set_num_threads(1)
    src, dst, kind, seed = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
    X = np.load(src)["X"]
    rng = np.random.default_rng(seed)
    torch.manual_seed(seed)
    n, p = X.shape
    df = pd.DataFrame(X, index=[f"c{i}" for i in range(n)], columns=[f"p{j}" for j in range(p)])
    df.index.name, df.columns.name = "Sample ID", "protein"
    long = df.stack().rename("intensity")
    val_idx = rng.random(len(long)) < 0.05
    train, val = long[~val_idx], long[val_idx]
    out = tempfile.mkdtemp()
    import contextlib
    import io

    with contextlib.redirect_stdout(io.StringIO()):
        if kind == "cf":
            from pimmslearn.sklearn.cf_transformer import CollaborativeFilteringTransformer

            m = CollaborativeFilteringTransformer(target_column="intensity", sample_column="Sample ID",
                                                  item_column="protein", n_factors=30, out_folder=out, batch_size=4096)
            m.fit(train, val, epochs_max=20, cuda=False, patience=1)
            full = m.transform(train)  # predicts every (cell, protein) absent from the observed long series
            Y = full.unstack().reindex(index=df.index, columns=df.columns).to_numpy(float)
        else:
            from pimmslearn.sklearn.ae_transformer import AETransformer

            tr = train.unstack().reindex(index=df.index, columns=df.columns)
            keep = tr.columns[tr.notna().sum(0) >= 2]  # PIMMS' scaler drops all-NaN columns; the rest stay NaN
            tr = tr[keep]
            va = val.unstack().reindex(index=df.index, columns=keep)
            m = AETransformer(hidden_layers=[64], latent_dim=10, out_folder=out, model="VAE", batch_size=10)
            m.fit(tr, va, epochs_max=300, cuda=False, patience=50)
            Y = m.transform(tr).reindex(index=df.index, columns=df.columns).to_numpy(float)
    np.savez_compressed(dst, Y=Y)


if __name__ == "__main__":
    main()
