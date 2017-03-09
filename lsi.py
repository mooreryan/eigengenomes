import logging, sys
from gensim import corpora, models, similarities

logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s',
                    level=logging.INFO)

if len(sys.argv) < 4:
    print "USAGE: python %s <outdir> <kmer_tfidf_freq.mm> <*.hash_bucket_counts (used in previous step)>" % sys.argv[0]
    exit(1)

fname = sys.argv[2]
corpus = corpora.MmCorpus(fname)

outdir = sys.argv[1]

# get the names of files used in previous step
sample_map = {}
for idx, f in enumerate(sys.argv[3:]):
    sample_map[idx] = f


# num_topics=20 is taken from streaming_eigenhashes.py
lsi = models.LsiModel(corpus, num_topics=20, id2word=sample_map, distributed=False, chunksize=200000)
lsi.save(outdir + '/kmer_lsi.gensim')
