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
    }
}
