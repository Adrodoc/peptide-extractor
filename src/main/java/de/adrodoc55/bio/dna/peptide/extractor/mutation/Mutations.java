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

import static java.lang.Character.toUpperCase;

import de.adrodoc55.bio.dna.AminoAcid;
import de.adrodoc55.bio.dna.peptide.extractor.UnknownAminoAcidException;
import de.adrodoc55.bio.dna.peptide.extractor.ValidationException;

/**
 * @author Adrodoc55
 */
public class Mutations {
  public static Mutation parse(String header)
      throws UnknownAminoAcidException, ValidationException {
    Mutation termination = Termination.parse(header);
    if (termination != null) {
      // Ignore terminating mutations, because they dont cause an aminoacid sequence alternation
      return null;
    }
    Mutation deletion = Deletion.parse(header);
    if (deletion != null) {
      return deletion;
    }
    Mutation insertion = Insertion.parse(header);
    if (insertion != null) {
      return insertion;
    }
    Mutation snp = SingleNucleotidePolymorphism.parse(header);
    if (snp != null) {
      return snp;
    }
    return null;
  }

  /**
   *
   * @param protein
   * @param mutationIndex (1 based)
   * @param expected
   * @throws ValidationException
   */
  public static void checkAmino(CharSequence protein, int mutationIndex, AminoAcid expected)
      throws ValidationException {
    char actualAmino = toUpperCase(protein.charAt(mutationIndex - 1));
    char expectedAmino = expected.getCharCode();
    checkValid(actualAmino != expectedAmino, "Incorrect aminoacid at index " + mutationIndex
        + "! Expected " + expectedAmino + " but got " + actualAmino);
  }

  public static void checkValid(boolean b, String message) throws ValidationException {
    if (b) {
      throw new ValidationException(message);
    }
  }
}
