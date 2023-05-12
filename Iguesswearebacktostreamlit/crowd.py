from bioservices import biomart

# Connect to the biomart server
server = biomart.BioMart(host='http://www.ensembl.org/biomart')

# Choose the dataset and attributes
dataset = "hsapiens_snp"
attributes = ["refsnp_id", "snp"]

# Choose the filters
filters = {
    "snp_filter": "rs1121980",
    "upstream_flank": 10,
    "downstream_flank": 10
}

# Create a query object
query = server.query(dataset=dataset, attributes=attributes, filters=filters)

# Execute the query and retrieve the results
results = query.results()

# Print the results
for result in results:
    print(result["refsnp_id"], result["snp"])