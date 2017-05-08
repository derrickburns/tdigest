# tdigest

This is an implementation of Ted Dunning's [Merging T-Digest](https://github.com/tdunning/t-digest/).

This implementation batches all inserts, including merges.  This is an improvement over the original implementation where merges require sorting of all points on each merge. 

This implementation does not support storing the incoming data with each centroid, (since obviously that is for testing).


