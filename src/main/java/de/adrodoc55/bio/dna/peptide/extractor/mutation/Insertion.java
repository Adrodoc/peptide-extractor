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

import static de.adrodoc55.bio.dna.peptide.extractor.mutation.Mutations.checkAmino;
import static de.adrodoc55.bio.dna.peptide.extractor.mutation.Mutations.checkValid;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.google.common.collect.ImmutableList;

import de.adrodoc55.bio.dna.AminoAcid;
import de.adrodoc55.bio.dna.peptide.extractor.UnknownAminoAcidException;
import de.adrodoc55.bio.dna.peptide.extractor.ValidationException;

/**
 * @author Adrodoc55
 */
public class Insertion implements Mutation {
  private static Pattern PATTERN =
      Pattern.compile(">.*:p\\.([A-Za-z]{3})(\\d+)_([A-Za-z]{3})(\\d+)ins((?:[A-Za-z]{3})+)$");

  public static Insertion parse(CharSequence header)
      throws ValidationException, UnknownAminoAcidException {
    Matcher matcher = PATTERN.matcher(header);
    if (matcher.find()) {
      AminoAcid aminoBeforeInsertion = AminoAcid.from3LetterCode(matcher.group(1));
      int indexBeforeInsertion = Integer.parseInt(matcher.group(2));
      AminoAcid aminoAfterInsertion = AminoAcid.from3LetterCode(matcher.group(3));
      int indexAfterInsertion = Integer.parseInt(matcher.group(4));

      checkValid(indexAfterInsertion == indexBeforeInsertion + 1,
          "Expected index after insertion (" + indexAfterInsertion
              + ") to be 1 greater than index before insertion (" + indexBeforeInsertion + ")");

      String insertedAminosString = matcher.group(5);
      List<AminoAcid> insertedAminos = new ArrayList<>(insertedAminosString.length() / 3);
      for (int i = 0; i + 2 < insertedAminosString.length(); i += 3) {
        AminoAcid insertedAmino =
            AminoAcid.from3LetterCode(insertedAminosString.substring(i, i + 3));
        insertedAminos.add(insertedAmino);
      }

      return new Insertion(indexBeforeInsertion, aminoBeforeInsertion, aminoAfterInsertion,
          insertedAminos);
    } else {
      return null;
    }
  }

  private final int indexBeforeInsertion;
  private final AminoAcid aminoBeforeInsertion;
  private final AminoAcid aminoAfterInsertion;
  private final ImmutableList<AminoAcid> insertedAminos;

  private Insertion(int indexBeforeInsertion, AminoAcid aminoBeforeInsertion,
      AminoAcid aminoAfterInsertion, Iterable<? extends AminoAcid> insertedAminos) {
    this.indexBeforeInsertion = indexBeforeInsertion;
    this.aminoBeforeInsertion = aminoBeforeInsertion;
    this.aminoAfterInsertion = aminoAfterInsertion;
    this.insertedAminos = ImmutableList.copyOf(insertedAminos);
  }

  public int getIndexBeforeInsertion() {
    return indexBeforeInsertion;
  }

  public AminoAcid getAminoBeforeInsertion() {
    return aminoBeforeInsertion;
  }

  public int getIndexAfterInsertion() {
    return indexBeforeInsertion + 1 + insertedAminos.size();
  }

  public AminoAcid getAminoAfterInsertion() {
    return aminoAfterInsertion;
  }

  public ImmutableList<AminoAcid> getInsertionAminos() {
    return insertedAminos;
  }

  @Override
  public CharSequence extractFromProtein(CharSequence protein, int offset)
      throws ValidationException {
    validateProtein(protein);
    int beginIndex = Math.max(0, getIndexBeforeInsertion() - offset);
    int endIndex = Math.min(protein.length(), getIndexAfterInsertion() - 1 + offset);
    return protein.subSequence(beginIndex, endIndex);
  }

  @Override
  public void validateProtein(CharSequence protein) throws ValidationException {
    checkAmino(protein, indexBeforeInsertion, aminoBeforeInsertion);
    int insertionSize = insertedAminos.size();
    int firstMutationIndex = indexBeforeInsertion + 1;
    for (int i = 0; i < insertionSize; i++) {
      checkAmino(protein, firstMutationIndex + i, insertedAminos.get(i));
    }
    checkAmino(protein, getIndexAfterInsertion(), aminoAfterInsertion);
  }

  @Override
  public String getUniqueSolution(CharSequence output) {
    return "INS-" + output;
  }
}
