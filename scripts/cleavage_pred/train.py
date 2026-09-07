from protai.annotation.data import DatasetGenerator
import pickle
import json
import os
import re
import matplotlib.pyplot as plt


def train(data_path, json_path):
	# Train split percent
	# Training data path
	# Save directory
	# Batch size
	dg = DatasetGenerator(0.9, json_path, data_path, bs=128)
	dg.run(100) # Number of epochs


def test(data_path, raw_data_path):
	#data_path = Path(data_path)
	model_path = data_path + "/model.p"
	vocab_path = data_path + "/vocab.pkl"
	datasplit_path = data_path + "/datasplit.json"

	results = DatasetGenerator.evaluate_later(model_path, vocab_path, datasplit_path, raw_data_path)
	outpath = data_path + "/tested.json"
	with open(outpath,"w") as fp:
		json.dump(results, fp, indent=4)

def predict(data_path, sequences):
	#data_path = Path(data_path)
	model_path = data_path + "/model.p"
	vocab_path = data_path + "/vocab.pkl"

def evaluate_test(test_path, outpath):
	with open(test_path, 'r') as fin:
		l_test = json.load(fin)
	l_diff = []
	for bac in l_test:
		name = bac['name'].lstrip('start-').rstrip('-stop')
		sequence = bac['sequence'].lstrip('start-').rstrip('-stop')
		labels = bac['labels'].lstrip('start-').rstrip('-stop')
		prediction = bac['prediction'].lstrip('start-').rstrip('-stop')
		labels = re.sub(r'before', 'b', labels)
		labels = re.sub(r'prop', 'p', labels)
		labels = re.sub(r'-', '', labels)
		prediction = re.sub(r'before', 'b', prediction)
		prediction = re.sub(r'prop', 'p', prediction)
		prediction = re.sub(r'-', '', prediction)
		label_leader = len(re.search('(^b.+?)p', labels).group(1))
		pred_leader = len(re.search('(^b.+?)p', prediction).group(1))
		difference = pred_leader - label_leader
		l_diff.append(difference)

	# Count the frequency of each integer
	unique_integers = sorted(set(l_diff))
	frequency = [l_diff.count(x) for x in unique_integers]
	# plot barplot
	for i, freq in enumerate(frequency):
		plt.text(unique_integers[i], freq, str(freq), ha='center', va='bottom')  
	plt.bar(unique_integers, frequency)
	plt.xlabel('Difference in cleavage site prediction')
	plt.ylabel('Counts')
	# plt.title('Bar Plot of Integers')
	plt.grid(False)  # Remove the grid
	plt.show()
	plt.savefig(os.path.join(outpath, 'prediction.test.png'))
	plt.savefig(os.path.join(outpath, 'prediction.test.pdf'))



if __name__ == "__main__":
	current_path = os.path.split(os.path.realpath(__file__))[0]
	data_path = os.path.join(current_path, 'training_data/annotation')
	json_path = data_path + "/train_set.json"
	test_path = data_path + "/tested.json"

	train(data_path, json_path)
	test(data_path, json_path)
	evaluate_test(test_path, data_path)
