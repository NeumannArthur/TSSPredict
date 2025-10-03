import pandas as pd
import torch
from evo2 import Evo2
import pandas as pd


# def load_and_label(csv_path, label, strand):
#     df = pd.read_csv(csv_path, sep="\t")
#     df["label"] = label
#     return df

# # Very short script just to combine Control strain as exemplary dataset -> not permanent solution!!!
# plus_tss = load_and_label("/ceph/ibmi/it/thesis_data/neumannarthur/results/Control/plus/[Window.Control]_plusStrand.tsv", strand="+", label=1)
# plus_neg = load_and_label("/ceph/ibmi/it/thesis_data/neumannarthur/results/Control/plus/[Window.Control]_plusStrand_negative.tsv", strand="-", label=0)
# minus_tss = load_and_label("/ceph/ibmi/it/thesis_data/neumannarthur/results/Control/minus/[Window.Control]_minusStrand.tsv", strand="+", label = 1)
# minus_neg = load_and_label("/ceph/ibmi/it/thesis_data/neumannarthur/results/Control/minus/[Window.Control]_minusStrand_negative.tsv", strand="-", label=0)

# master_df = pd.concat([plus_tss, minus_tss, plus_neg, minus_neg])
# # master_df = master_df.sample(frac=1, random_state=42).reset_index(drop=True)

# print(master_df.head())
# master_df.to_csv("/ceph/ibmi/it/thesis_data/neumannarthur/figures/joined_db.tsv", sep="\t")

def create_embeddings(layer_num=26):
    model = Evo2('evo2_7b_base')

    positions = pd.read_csv('/ceph/ibmi/it/thesis_data/neumannarthur/figures/joined_db.tsv', sep='\t')[["Bases", "contigPos"]]

    def get_window_embedding(sequence, model, layer_name=f'blocks.{layer_num}'):
        input_ids = torch.tensor(
            model.tokenizer.tokenize(sequence),
            dtype=torch.int,
        ).unsqueeze(0).to('cuda:0')
        with torch.no_grad():
            _, embeddings = model(input_ids, return_embeddings=True, layer_names=[layer_name])
        return embeddings[layer_name][0, -1, :].cpu().to(torch.float32).numpy()


    positions["Embedding"] = positions["Bases"].apply(
        lambda seq: get_window_embedding(seq, model).tolist()
    )

    pd.merge(pd.read_csv('/ceph/ibmi/it/thesis_data/neumannarthur/figures/joined_db.tsv', sep='\t'), 
            positions[["contigPos", "Embedding"]], 
            on="contigPos", 
            how="left").to_csv(f'/ceph/ibmi/it/thesis_data/neumannarthur/data/ecoli_data/features/[Embeddings.Control]layer_{layer_num}.tsv')

    print(pd.read_csv(f'/ceph/ibmi/it/thesis_data/neumannarthur/data/ecoli_data/features/[Embeddings.Control]layer_{layer_num}.tsv').head(50))

for i in range(22, 33, 2):
    create_embeddings(i)