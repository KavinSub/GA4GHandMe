import json
import requests
from ga4gh.client import client
from collections import defaultdict

c = client.HttpClient("http://1kgenomes.ga4gh.org")
dataset = c.search_datasets().next()

for variant_set in c.search_variant_sets(dataset_id=dataset.id):
	if variant_set.name == "phase3-release":
		variant_set_id = variant_set.id

# Fetch all call_set_ids
call_set_ids = []
for call_set in c.search_call_sets(variant_set_id=variant_set_id):
	call_set_ids.append(call_set.id)

def is_equal_variant(x, y):
	return x["reference_bases"] == y[0] and x["alternate_bases"][0] == y[1]

def nearest_relative(data):
	print("####")
	counter = defaultdict(int)
	for index, variant in enumerate(data):
		print(variant)
		for variant_with_calls in c.search_variants(call_set_ids=call_set_ids, variant_set_id=variant_set_id, reference_name="1", start=variant["start"], end=variant["end"]):
			ref_base = variant_with_calls.reference_bases
			alt_base = str(variant_with_calls.alternate_bases[0])
			if is_equal_variant(variant, (ref_base, alt_base)):
				for call in variant_with_calls.calls:
					if call.genotype[0] | call.genotype[1]:
						counter[call.call_set_id] += 1
	nearest_relative = max(counter.keys(), key=lambda x: counter[x])
	print("Call set id of nearest relative: {}\nVariations in common: {}".format(nearest_relative, counter[nearest_relative]))
	print("####")