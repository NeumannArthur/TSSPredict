from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report
import pandas as pd
import numpy as np
import ast 
import sys

print("<------------------------------ Loading data... ------------------------------>")

data = pd.read_csv(sys.argv[1])

print("<------------------------------ Converting data... ------------------------------>")
data["coverage"] = data["coverage"].apply(lambda x: np.array(ast.literal_eval(x), dtype=float).ravel())
data["Embedding"] = data["Embedding"].apply(lambda x: np.array(ast.literal_eval(x), dtype=float).ravel())

print(data["Embedding"].head(50))

cov_df = pd.DataFrame(data["coverage"].tolist(), index=data.index)
emb_df = pd.DataFrame(data["Embedding"].tolist(), index=dsxata.index)

print("<------------------------------ Splitting data... ------------------------------>")

X_data = pd.concat([data[["contigPos", "windowBegin", "windowEnd"]], cov_df, emb_df], axis=1)
# X_data = X_data.apply(pd.to_numeric, errors="coerce")
y_data = data["label"].astype(int)

X_train, X_test, y_train, y_test = train_test_split(
    X_data, y_data, test_size=0.2, random_state=42, stratify=y_data
)


# # print(y_data.head(50))
# # print(X_data.head(50))

print("<------------------------------ Training model... ------------------------------>")

clf = RandomForestClassifier(n_estimators=200, random_state=42, n_jobs=-1)
clf.fit(X_train.to_numpy(), y_train.to_numpy())

y_pred = clf.predict(X_test.to_numpy())

print("Accuracy:", accuracy_score(y_test, y_pred))
print(classification_report(y_test, y_pred))

