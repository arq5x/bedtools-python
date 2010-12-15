#!/usr/bin/env python

import unittest
import os
from bedtools import IntervalFile
PATH = os.path.dirname(__file__)

class IntervalFileTest(unittest.TestCase):
    file = "data/rmsk.hg18.chr21.bed"
    def setUp(self):
        self.file = os.path.join(PATH, self.file)
        self.bed = IntervalFile(self.file)

    def testOverlaps(self):
        hits = self.bed.search("chr21", 9719768, 9739768)
        print len(hits)
        self.assertEqual(len(hits), 8)
        for hit in hits:
            self.assert_(hit.start <= 9739768 and hit.end >= 9719768)

    def testStrands(self):
        hits = self.bed.search("chr21", 9719768, 9739768, "+")
        for hit in hits:
            self.assert_(hit.strand == '+')
        hits = self.bed.search("chr21", 9719768, 9739768, "-")
        for hit in hits:
            self.assert_(hit.strand == '-')

class IntervalFileGzTest(IntervalFileTest):
      file = "data/rmsk.hg18.chr21.bed.2.gz"


if __name__ == "__main__":
    unittest.main()


