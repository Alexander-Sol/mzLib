using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public static class FastaScratchPad
    {

        [Test]
        public static void ArabidaDigest()
        {
            var proteins = ProteinDbLoader.LoadProteinFasta(@"C:\Users\Alex\Documents\Proteomes\Arabidopsis_Thaliana.fasta",
                true, DecoyType.Reverse, isContaminant: false, out var errors);

            List<string> peptideBaseSequences = new();

            foreach(Protein protein in proteins)
            {
                protein.Digest(new DigestionParams(), null, null);
                peptideBaseSequences.AddRange(protein.Digest(new DigestionParams(), null, null).Select(p => p.BaseSequence));
            }

            using (StreamWriter writer = new(@"C:\Users\Alex\Documents\Proteomes\Arabidopsis_Thaliana_tryptic_Peptides.tsv"))
            {
                writer.WriteLine("Peptide Base Sequence");
                foreach(string seq in peptideBaseSequences)
                {
                    writer.WriteLine(seq);
                }
            }
        }

        [Test]
        public static void HumanReferenceDigest()
        {
            var proteins = ProteinDbLoader.LoadProteinXML(@"C:\Users\Alex\Documents\Proteomes\Uniprot_H_Sapiens_Reviewed_11_28_23.xml",
                true, DecoyType.Reverse, null, isContaminant: false, null, out var unknownMods);

            List<PeptideWithSetModifications> peptides = new();

            foreach (Protein protein in proteins)
            {
                peptides.AddRange(protein.Digest(new DigestionParams(), null, null));
            }

            using (StreamWriter writer = new(@"C:\Users\Alex\Documents\Proteomes\H_Sapiens_tryptic_Peptides.tsv"))
            {
                writer.WriteLine("PeptideModSeq\tPeptideMonoisotopicMass");
                foreach (var peptide in peptides)
                {
                    if (peptide.BaseSequence.Contains('U')) continue;
                    writer.WriteLine(peptide.BaseSequence + '\t' + 
                        peptide.MonoisotopicMass.ToString());
                }
            }
        }
    }
}
