# from protai.classification.data import DatasetGenerator as CDG
from protai.annotation.data import DatasetGenerator as ADG
from pathlib import Path

#models_dir = Path("/home/ygao/models") # downloaded from releases! 
models_dir = Path("/home/dwzhang/test/gy/training_data/")


#class_model_dir = models_dir / "classification"
#class_model_path = class_model_dir / "model.p"
#class_vocab_path = class_model_dir / "vocab.pkl"
annot_model_dir = models_dir / "annotation"
annot_model_path = annot_model_dir / "model.p"
annot_vocab_path = annot_model_dir / "vocab.pkl"

sequences = [
    {
        "sequence": "MISSHQKTLTDKELALISGGKTHYPTNAWKSLWKGFWESLRYTDGF",
        "name": "unique_name1"},
    {
    "sequence": "MEEIISFEFSNGSKLLDSNVYEYKVGKEGTIKIDFPDSNGFVAVHKNNTETAYIKAPYMKFCTGRHA",
    "name": "GCF_000143745.1_+_NZ_GL379763.1_390"},
    {
    "sequence": "MISSHQKTLTDKELALISGGKTHYPTNAWKSLWKDDDESLRYTDGFTTTTTFFGSAAAAAAAAAAKKKKKKKKKKKKKKKKKKKKKKKTTTTTTTTTTTTTTTTTTTTTTMMMMMMMMMSSSSSSSSSSYYYYYYYYYYYYDDDDDDD",
    "name": "unique_name3"}
]

#class_predictions = CDG.predict(class_model_path, class_vocab_path, sequences)
cleavage_predictions = ADG.predict(annot_model_path, annot_vocab_path, sequences)
'''
import json
print("Class predictions")
print(json.dumps(class_predictions, indent=4))

print("Cleavage predictions")
print(json.dumps(cleavage_predictions, indent=4))
'''
print(cleavage_predictions)
