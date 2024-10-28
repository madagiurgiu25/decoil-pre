import unittest
from unittest.mock import patch

import numpy as np
import bionumpy
import vcfpy
from vcfpy import InvalidRecordException, Record, BreakEnd, Call

from decoil.encode.encode import parsevcf
from decoil.utils import VCF_PROP as vp
from decoil.encode.operations import pass_filter

class TestInput(unittest.TestCase):

    def test_header_sniffles2(self):
        """Test correct format for sniffles2"""
        try:
            svinfo, collection_breakpoints, _ = parsevcf("tests/examples/ecdna1/sniffles2.vcf",vp.SNIFFLES2)
        except Exception as e:
            self.fail(f"Unexpected exception raised: {e}")
        # Optionally, check that the result is correct
        # svinfo[chr2_15585355] [('Sniffles2.BND.2S1', 'chr3', '11049996', 'BND', 47, 88, '0/1', '+-'), ('Sniffles2.BND.3S1', 'chr3', '11056099', 'BND', 55, 80, '0/1', '-+')]
        self.assertEqual(svinfo["chr2_15585355"][0][2], '11049996')
        self.assertEqual(svinfo["chr2_15585355"][1][1], 'chr3')
        self.assertEqual(svinfo["chr2_15585355"][1][7], '-+')
        
         # collection_breakpoints: defaultdict(<class 'list'>, {'chr2': [15585355, 15585355, 15633375, 16521051, 15585355], 'chr3': [11049996, 11056099, 10981201], 'chr12': [68807720, 68970909]})
        self.assertEqual(collection_breakpoints["chr2"], [15585355, 15585355, 15633375, 16521051, 15585355])
        
    def test_header_sniffles1(self):
        """Test correct format for sniffles1"""
        try:
            svinfo, collection_breakpoints, _ = parsevcf("tests/examples/ecdna1/sniffles1.vcf",vp.SNIFFLES)
            print(svinfo)
            print(collection_breakpoints)
        except Exception as e:
            self.fail(f"Unexpected exception raised: {e}")
        
        # 'chr2_15585355': [('4', 'chr3', '10981201', 'TRA', 24, 0, '1/1', '--'), ('5', 'chr3', '11060000', 'TRA', 123, 0, '1/1', '-+')]
        self.assertEqual(svinfo["chr2_15585355"][0][2], '10981201')
        self.assertEqual(svinfo["chr2_15585355"][1][1], 'chr3')
        self.assertEqual(svinfo["chr2_15585355"][1][7], '-+')
        
        #collection_breakpoints["chr2"]: [15633375, 16521051, 15585355, 15585355, 15633377, 16628304]  
        self.assertEqual(collection_breakpoints["chr2"], [15633375, 16521051, 15585355, 15585355, 15633377, 16628304])
    
    def test_header_cutesv(self):
        """Test correct format for cutesv"""
        try:
            svinfo, collection_breakpoints, _ = parsevcf("tests/examples/ecdna1/cutesv.vcf",vp.CUTESV)
            print(svinfo)
            print(collection_breakpoints)
        except Exception as e:
            self.fail(f"Unexpected exception raised: {e}")
        
        # 'chr2_15585356':[('cuteSV.BND.0', 'chr3', '10981202', 'BND', 24, 0, '1/1', '++')]
        self.assertEqual(svinfo["chr2_15585356"][0][2], '10981202')
        self.assertEqual(svinfo["chr2_15585356"][0][1], 'chr3')
        self.assertEqual(svinfo["chr2_15585356"][0][7], '++')
        
        #collection_breakpoints["chr2"]: [15585356, 15585357, 15585359, 15633375, 15633376, 16628305]  
        self.assertEqual(collection_breakpoints["chr2"], [15585356, 15585357, 15585359, 15633375, 15633376, 16628305])
    
    def test_header_nanomonsv(self):
        """Test correct format for nanomonsv"""
        try:
            svinfo, collection_breakpoints, _ = parsevcf("tests/examples/ecdna1/cutesv.vcf",vp.NANOMONSV)
            print(svinfo)
            print(collection_breakpoints)
        except Exception as e:
            self.fail(f"Unexpected exception raised: {e}")
        
        # 'chr2_15585356':[('cuteSV.BND.0', 'chr3', '10981202', 'BND', 24, 0, '1/1', '++')]
        self.assertEqual(svinfo["chr2_15585356"][0][2], '10981202')
        self.assertEqual(svinfo["chr2_15585356"][0][1], 'chr3')
        self.assertEqual(svinfo["chr2_15585356"][0][7], '++')
        
        #collection_breakpoints["chr2"]: [15585356, 15585357, 15585359, 15633375, 15633376, 16628305]  
        self.assertEqual(collection_breakpoints["chr2"], [15585356, 15585357, 15585359, 15633375, 15633376, 16628305])
    
    @patch('sys.exit')  # Mock sys.exit so it doesn't actually exit the interpreter
    def test_header_sniffles2_wrong_with_sniffles1(self, mock_exit):
        """Test wrong assignment of file format"""
        
        # This triggers an error
        parsevcf("tests/examples/ecdna1/sniffles1.vcf","sniffles2")
        mock_exit.assert_called_with(1)
        
    @patch('sys.exit')  # Mock sys.exit so it doesn't actually exit the interpreter
    def test_header_sniffles2_wrong_with_cutesv(self, mock_exit):
        """Test wrong assignment of file format"""
        
        # This triggers an error
        parsevcf("tests/examples/ecdna1/cutesv.vcf","sniffles1")
        mock_exit.assert_called_with(1)
        
    def test_none_DR(self):
        with self.assertRaises(Exception) as context:
            record  = Record('chr2', 15585355, ['Sniffles2.BND.1S1'], 'N', [BreakEnd('chr3', 10981201, '+', '-', 'N', True)], 55, ['GT'], {'PRECISE': True, 'SVTYPE': 'BND', 'SUPPORT': 14, 'COVERAGE': [0.0, 0.0, 136.0, 134.0, 132.0], 'STRAND': '+-', 'AF': 0.104, 'CHR2': 'chr3', 'STDEV_POS': 0.463}, ['GT', 'GQ', 'DR', 'DV'], [Call('SAMPLE', {'GT': '0/0', 'GQ': 60, 'DR': None, 'DV': 14})])
            pass_filter(record)
        self.assertEqual(str(context.exception), "DR has not value. Your VCF might not be genotyped. Rerun SV calling using --genotype")
    
    def test_none_DV(self):
        with self.assertRaises(Exception) as context:
            record  = Record('chr2', 15585355, ['Sniffles2.BND.1S1'], 'N', [BreakEnd('chr3', 10981201, '+', '-', 'N', True)], 55, ['GT'], {'PRECISE': True, 'SVTYPE': 'BND', 'SUPPORT': 14, 'COVERAGE': [0.0, 0.0, 136.0, 134.0, 132.0], 'STRAND': '+-', 'AF': 0.104, 'CHR2': 'chr3', 'STDEV_POS': 0.463}, ['GT', 'GQ', 'DR', 'DV'], [Call('SAMPLE', {'GT': '0/0', 'GQ': 60, 'DR': 121, 'DV': None})])
            pass_filter(record)
        self.assertEqual(str(context.exception), "DV has not value. Your VCF might not be genotyped. Rerun SV calling using --genotype")

    # def test_vcf(self):
    #     file1 = "tests/examples/vcfs/sim1010_sniffles1.vcf"
    #     file2 = "tests/examples/vcfs/sim1010_sniffles2.vcf"
        
    #     print("pyvcf and sniffles2----")
    #     # read using vcf reader
    #     reader = vcfpy.Reader.from_path(file2)
    #     for record in reader:
    #         print(record)
        
    #     # print()
    #     # print("pyvcf and sniffles 2-------",flush=True)
        
    #     # reader = vcfpy.Reader.from_path(file2)
    #     # for record in reader:
    #     #     if record.INFO[vp.SVTYPE] == 'BND':
    #     #         print(record)
        
    #     print()
    #     print("bionumpy and sniffles2----")
    #     f = bionumpy.open(file2)
    #     for chunk in f.read_chunks():
    #         # we can iterate over the entries
    #         for single_entry in chunk.to_iter():
    #             print(single_entry)
                
        
            
    #     self.assertEqual("True","False")
        

if __name__ == '__main__':
	unittest.main()
