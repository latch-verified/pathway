from Bio import Entrez


with Entrez.efetch(db="gene", id="51146", retmode="xml") as handle:
    print(handle)
    # record = Entrez.read(handle)

# print(record)
