from collections import defaultdict
from gensim import corpora, models, similarities
from scipy.spatial import distance

import logging, sys
import numpy as np
import random
import os

# TODO clean up logging

logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s',
                    level=logging.INFO)

if len(sys.argv) != 3:
    print "USAGE: python %s <outdir> <kmer_tfidf_freq.mm>" % sys.argv[0]
    exit(1)


# To quickly adjust params for testing.
try:
    large = os.environ["LARGE"]
except KeyError:
    large = False

if large == "t":
    # TODO depending on the size of the corpus, you may not generate any chunks!


    print "It's a LARGE doc"
    # TODO: This should be a tunable param
    CLUSTER_THRESHOLD = 0.8

    # this is ~0.004 from the paper
    PERCENTAGE_OF_KMER_DOCS_FOR_SEEDING = 0.004

    # this is ~10,000 from the paper
    KMER_DOC_CHUNK_SIZE = 10000
else:
    print "It's a SMALL doc"
    # For the small test dataset
    CLUSTER_THRESHOLD = 0.8
    PERCENTAGE_OF_KMER_DOCS_FOR_SEEDING = 0.5
    KMER_DOC_CHUNK_SIZE = 5



outdir = sys.argv[1]
fname = sys.argv[2]
corpus = corpora.MmCorpus(fname)

# get the names of files used in previous step, for id2word param (not actually needed)
# sample_map = {}
# for idx, f in enumerate(sys.argv[3:]):
#     sample_map[idx] = f


# # A magic number from the original program
num_topics = round(len(sys.argv[3:]) * 0.8)

# set it to at least 5 (for testing)
if num_topics < 5:
    num_topics = 5

# from paper
READING_CHUNK_SIZE = 200000

lsi = models.LsiModel(corpus,
                      num_topics=num_topics,
                      #id2word=sample_map,
                      distributed=False,
                      chunksize=READING_CHUNK_SIZE)
# lsi.save(outdir + '/kmer_lsi.gensim')

# Done training the lsi now we cluster kmers

# map the kmer_docs into topic space
corpus_transform = lsi[corpus]


num_kmer_docs = len(corpus_transform)
num_kmer_docs_to_sample = num_kmer_docs * PERCENTAGE_OF_KMER_DOCS_FOR_SEEDING

NUM_CHUNKS = int(round((num_kmer_docs * PERCENTAGE_OF_KMER_DOCS_FOR_SEEDING) / KMER_DOC_CHUNK_SIZE))

kmer_doc_seeds = np.zeros((0, num_topics))
cluster_sizes = {}
kmer_doc_clusters = defaultdict(list)

kmer_docs = []
total_sampled = 0
# TODO our random kmer chunks will be near one another...is you wont
# likely get a kmer from the beginning and from the end of the
# corpus...does it matter?
i = 0
chunkno = 0
for kmer_doc in corpus_transform:
    i += 1
    if (i % 100000) == 0:
        print "processing doc: %s of %s" % (i, num_kmer_docs)

    if chunkno > NUM_CHUNKS:
        break

    if len(kmer_docs) < KMER_DOC_CHUNK_SIZE:
        # should we sample this one?
        if random.random() < PERCENTAGE_OF_KMER_DOCS_FOR_SEEDING:
            total_sampled += 1
            kmer_docs.append(kmer_doc)
    else:
        chunkno += 1
        print "Processing Chunk: %s" % chunkno

        for kmer_doc in kmer_docs:
            # kmer_doc has form [(topic0, val), (topic1, val), ..., (topic num_topics-1, val)]
            topic_values = np.zeros(num_topics)

            # get the values for this kmer_doc into an ary
            for topic, val in kmer_doc:
                topic_values[topic] = val

            num_kmer_doc_seeds = kmer_doc_seeds.shape[0]
            if num_kmer_doc_seeds > 0:
                # find the distance between this kmer_doc's values and the kmer_doc seeds
                dist = distance.cdist([topic_values], kmer_doc_seeds, 'cosine')
                # the index (cluster number) of the smallest distance
                clust = dist.argsort()[0][0]

                # if it is close enough
                if dist[0][clust] < 1 - CLUSTER_THRESHOLD:
                    # add kmer_doc to this cluster
                    kmer_doc_clusters[clust].append(topic_values)
                    cluster_sizes[clust] += 1
                else:
                    # add a new seed
                    kmer_doc_seeds = np.concatenate((kmer_doc_seeds, [topic_values]))

                    idx = len(kmer_doc_clusters)
                    # make a new cluster, and add this kmer_doc to it
                    kmer_doc_clusters[idx].append(topic_values)

                    # this size is just a place holder, will increment later
                    cluster_sizes[idx] = 1
            else:
                # add it to the seeds
                kmer_doc_seeds = np.array([topic_values])

                idx = 0 # the first one

                # this is the first topic_values we've seen, so it becomes a cluster
                kmer_doc_clusters[idx].append(topic_values)

                # this size is just a place holder, will increment later
                cluster_sizes[idx] = 1

        for cluster, topic_values in kmer_doc_clusters.items():
            # get the mean values for all kmer_docs in this cluster
            mean_topic_values = np.array(topic_values).sum(0) / len(topic_values)

            # replace the old seed with the mean of the cluster (sort of a "consensus" seed of the cluster")
            kmer_doc_seeds[cluster, :] = mean_topic_values

            # increment the cluster size (TODO...is this correct?)
            # cluster_sizes[cluster] += len(topic_values)

        # only send the seeds to the next step as members of the
        # kmer_clusters. TODO -- this will save memory, but is it bad to
        # do? I think the orig does this as well.
        # TODO check this
        # TODO .keys() creates a list of all keys...could use too much mem
        for key in kmer_doc_clusters.keys():
            kmer_doc_clusters[key] = [kmer_doc_seeds[key]]

        # merge any clusters whose seed is within the cluster threshold
        seed_dist = distance.cdist(kmer_doc_seeds, kmer_doc_seeds, 'cosine')
        nrows = seed_dist.shape[0]
        ncols = seed_dist.shape[1]

        if nrows != ncols:
            print "ERROR"
            exit(111)

        # TODO would merging the clusters (like the paper says but the
        # orig program doesn't do), be better than removing the one of the
        # clusters?
        clusters_to_del = {}
        for row in range(nrows):
            for col in range(row+1, ncols):
                if seed_dist[row][col] < 1 - CLUSTER_THRESHOLD:
                    clusters_to_del[col] = row


        # for testing always delete col 2 and col 3
        # for the test Im assuming that (1, 2) and (2, 3) are both within the threshold.
        # SO this would want you to delete both 2 and 3, but then we'd lose a whole cluster....
        # if len(cluster_sizes) >= 4:
        #     print "IN THE TEST"
        #     print kmer_doc_seeds
        #     print kmer_doc_clusters
        #     print cluster_sizes
        #     clusters_to_del[2] = 1
        #     clusters_to_del[3] = 2
        #     print clusters_to_del

        cluster_sizes_before_merging = sum(cluster_sizes.values())
        new_kmer_doc_seeds = []
        new_kmer_doc_clusters = defaultdict(list)
        new_cluster_sizes = {}
        already_deleted = {}
        if len(clusters_to_del) > 0:
            # print "We will merge some clusters now"
            for cluster_idx in xrange(len(kmer_doc_clusters)):
                # clusters_to_del[cluster_idx] is the cluster within
                # threshold of the one to delete
                if (cluster_idx in clusters_to_del) and (clusters_to_del[cluster_idx] not in already_deleted):
                    # print "Merging %s into %s" % (cluster_idx, clusters_to_del[cluster_idx])
                    # delete this cluster
                    already_deleted[cluster_idx] = True

                    idx_to_keep = clusters_to_del[cluster_idx]
                    if idx_to_keep in new_cluster_sizes:
                        # add the count of the cluster that will be
                        # deleted to the one that will be kept
                        new_cluster_sizes[idx_to_keep] += cluster_sizes[cluster_idx]
                    else:
                        new_cluster_sizes[idx_to_keep] = cluster_sizes[cluster_idx]

                else:
                    # print "Keeping %s" % cluster_idx
                    new_idx = len(new_kmer_doc_clusters)

                    # TODO better might be to take the mean of the two
                    # that will be merged rather than just deleting
                    # one
                    new_kmer_doc_clusters[new_idx] = [kmer_doc_clusters[cluster_idx]]
                    new_kmer_doc_seeds.append(kmer_doc_seeds[cluster_idx])
                    new_cluster_sizes[new_idx] = cluster_sizes[cluster_idx]

            kmer_doc_clusters = new_kmer_doc_clusters
            kmer_doc_seeds = np.array(new_kmer_doc_seeds)
            cluster_sizes = new_cluster_sizes

            clusters_to_del = {}
            already_deleted = {}

            cluster_sizes_after_merging = sum(cluster_sizes.values())
            # an invariant
            if cluster_sizes_before_merging != cluster_sizes_after_merging:
                print "WARNING -- something wonky about the counts"
                # print cluster_sizes_before_merging
                # print cluster_sizes_after_merging
                # exit(112)
                # TODO fix this


        # reset it
        kmer_docs = []

# Write the cluster sizes
size_fname = fname + ".kmer_doc_cluster_sizes.txt"
size_f = open(size_fname, "w")
num_clusters = len(cluster_sizes)
size_f.write("%s\n" % num_clusters)
for cluster_idx, size in cluster_sizes.items():
    string = "%s %s\n" % (cluster_idx, size)
    size_f.write(string)

size_f.close()
print "Wrote %s" % size_fname

# Write the cluster seeds
seeds_fname = fname + ".kmer_doc_seeds.txt"
seeds_f = open(seeds_fname, "w")
nrows = len(kmer_doc_seeds)
ncols = len(kmer_doc_seeds[0])
seeds_f.write("%s %s\n" % (nrows, ncols))
for vals in kmer_doc_seeds:
    string = " ".join(str(val) for val in vals) + "\n"
    seeds_f.write(string)

seeds_f.close()
print "Wrote %s" % seeds_fname

# save it and serialize so we can have random access to the corpus
# TODO is there another way to do this? Cos it is slow
print "Serializing transformed corpus...this could take a while"
# This is needed for random access in the next step
transformed_corpus_fname = fname + ".lsi_transformed.mm"
corpora.MmCorpus.serialize(transformed_corpus_fname, corpus_transform)
print "Wrote %s" % transformed_corpus_fname

print "DONE!"
