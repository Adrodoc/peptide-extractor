/*
 * Peptide Extractor: A command line tool for transforming the output of the Ensembl Variant Effect
 * Predictor ProteinSeqs Plugin into NetMHC-readable peptide fragments that are affected by
 * mutation.
 *
 * © Copyright (C) 2017 Adrodoc55
 *
 * This file is part of Peptide Extractor.
 *
 * Peptide Extractor is free software: you can redistribute it and/or modify it under the terms of
 * the GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * Peptide Extractor is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with Peptide Extractor.
 * If not, see <http://www.gnu.org/licenses/>.
 *
 *
 *
 * Peptide Extractor: Ein Kommandozeilen Werkzeug um den Output des Ensembl Variant Effect Predictor
 * ProteinSeqs Plugins in NetMHC-lesbare Peptid Fragmente die von einer Mutation beeinflusst sind
 * umzuwandeln.
 *
 * © Copyright (C) 2017 Adrodoc55
 *
 * Diese Datei ist Teil von Peptide Extractor.
 *
 * Peptide Extractor ist freie Software: Sie können diese unter den Bedingungen der GNU General
 * Public License, wie von der Free Software Foundation, Version 3 der Lizenz oder (nach Ihrer Wahl)
 * jeder späteren veröffentlichten Version, weiterverbreiten und/oder modifizieren.
 *
 * Peptide Extractor wird in der Hoffnung, dass es nützlich sein wird, aber OHNE JEDE
 * GEWÄHRLEISTUNG, bereitgestellt; sogar ohne die implizite Gewährleistung der MARKTFÄHIGKEIT oder
 * EIGNUNG FÜR EINEN BESTIMMTEN ZWECK. Siehe die GNU General Public License für weitere Details.
 *
 * Sie sollten eine Kopie der GNU General Public License zusammen mit Peptide Extractor erhalten
 * haben. Wenn nicht, siehe <http://www.gnu.org/licenses/>.
 */
package de.adrodoc55.bio.dna.peptide.extractor.mutation;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import de.adrodoc55.bio.dna.Nucleotide;
import de.adrodoc55.bio.dna.peptide.extractor.UnknownNucleotideException;
import de.adrodoc55.bio.dna.peptide.extractor.ValidationException;

/**
 * @author Adrodoc55
 */
public class SilentSingleNucleotidePolymorphism implements Mutation {
  private static Pattern PATTERN =
      Pattern.compile(">.*:c\\.(\\d+)([A-Za-z])>([A-Za-z])\\(p\\.=\\)$");

  public static SilentSingleNucleotidePolymorphism parse(CharSequence header)
      throws UnknownNucleotideException {
    Matcher matcher = PATTERN.matcher(header);
    if (matcher.find()) {
      int mutationIndex = Integer.parseInt(matcher.group(1));
      Nucleotide nativeNucleotide = Nucleotide.fromCharCode(matcher.group(2).charAt(0));
      Nucleotide mutatedNucleotide = Nucleotide.fromCharCode(matcher.group(3).charAt(0));
      return new SilentSingleNucleotidePolymorphism(mutationIndex, nativeNucleotide,
          mutatedNucleotide);
    } else {
      return null;
    }
  }

  private final int mutationIndex;
  private final Nucleotide nativeNucleotide;
  private final Nucleotide mutatedNucleotide;

  public SilentSingleNucleotidePolymorphism(int mutationIndex, Nucleotide nativeNucleotide,
      Nucleotide mutatedNucleotide) {
    this.mutationIndex = mutationIndex;
    this.nativeNucleotide = nativeNucleotide;
    this.mutatedNucleotide = mutatedNucleotide;
  }

  public int getMutationIndex() {
    return mutationIndex;
  }

  public Nucleotide getNaitveNucleotide() {
    return nativeNucleotide;
  }

  public Nucleotide getMutatedNucleotide() {
    return mutatedNucleotide;
  }

  @Override
  public CharSequence extractFromProtein(CharSequence protein, int offset)
      throws ValidationException {
    throw new UnsupportedOperationException();
  }

  @Override
  public void validateProtein(CharSequence protein) throws ValidationException {}

  @Override
  public String getUniqueSolution(CharSequence output) {
    throw new UnsupportedOperationException();
  }
}
