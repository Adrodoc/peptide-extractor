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
package de.adrodoc55.bio.dna.peptide.extractor.main;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Set;
import java.util.regex.Matcher;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParameterException;
import com.google.common.base.Charsets;
import com.google.common.base.Joiner;
import com.google.common.collect.BiMap;
import com.google.common.collect.HashBiMap;
import com.google.common.collect.ImmutableList;
import com.google.common.io.Files;
import com.google.common.io.Resources;

/**
 * @author Adrodoc55
 */
public class PeptideExtractorMain {
  public static void main(String[] args) throws IOException {
    PeptideExtractorParameter params = new PeptideExtractorParameter();
    JCommander jc = new JCommander(params);
    jc.setProgramName("java -jar peptide-extractor.jar");
    try {
      jc.parse(args);
      if (params.isHelp()) {
        jc.usage();
        return;
      }
      URL inputUrl = params.getInput().toURI().toURL();
      ImmutableList<String> lines = Resources.asCharSource(inputUrl, Charsets.UTF_8).readLines();

      File outputFile = params.getOutput();
      Files.createParentDirs(outputFile);
      outputFile.createNewFile();

      List<String> outputLines = transform(lines, params);
      String output = Joiner.on("\n").join(outputLines);
      Files.write(output, outputFile, Charsets.UTF_8);
    } catch (ParameterException ex) {
      System.err.println(ex.getLocalizedMessage());
      System.err.println("Run with '-h' to print help");
    }
  }

  private static final BiMap<String, Character> AMINO_ACIDS = HashBiMap.create(20);
  static {
    AMINO_ACIDS.put("ALA", 'A');
    AMINO_ACIDS.put("ARG", 'R');
    AMINO_ACIDS.put("ASN", 'N');
    AMINO_ACIDS.put("ASP", 'D');
    AMINO_ACIDS.put("CYS", 'C');
    AMINO_ACIDS.put("GLN", 'Q');
    AMINO_ACIDS.put("GLU", 'E');
    AMINO_ACIDS.put("GLY", 'G');
    AMINO_ACIDS.put("HIS", 'H');
    AMINO_ACIDS.put("ILE", 'I');
    AMINO_ACIDS.put("LEU", 'L');
    AMINO_ACIDS.put("LYS", 'K');
    AMINO_ACIDS.put("MET", 'M');
    AMINO_ACIDS.put("PHE", 'F');
    AMINO_ACIDS.put("PRO", 'P');
    AMINO_ACIDS.put("SER", 'S');
    AMINO_ACIDS.put("THR", 'T');
    AMINO_ACIDS.put("TRP", 'W');
    AMINO_ACIDS.put("TYR", 'Y');
    AMINO_ACIDS.put("VAL", 'V');
  }

  private static List<String> transform(List<String> lines, PeptideExtractorParameter params) {
    List<String> result = new ArrayList<>();
    Set<String> uniqueSolutions = new HashSet<>();
    for (int i = 0; i < lines.size();) {
      int headerLine = i++;
      String header = lines.get(headerLine);
      if (header.startsWith(params.getHeaderPrefix())) {

        StringBuilder sb = new StringBuilder();
        while (i < lines.size() && !lines.get(i).startsWith(params.getHeaderPrefix())) {
          sb.append(lines.get(i++));
        }

        if (params.getHeaderPattern().matcher(header).matches()) {
          Matcher mutationIndexMatcher = params.getMutationIndexPattern().matcher(header);
          mutationIndexMatcher.find();
          // 1 based index
          int mutationIndex = Integer.parseInt(mutationIndexMatcher.group(1));

          Matcher nativeAminoMatcher = params.getNativeAminoPattern().matcher(header);
          nativeAminoMatcher.find();
          String nativeAmino = nativeAminoMatcher.group(1).toUpperCase(Locale.ENGLISH);

          Matcher mutatedAminoMatcher = params.getMutatedAminoPattern().matcher(header);
          mutatedAminoMatcher.find();
          String mutatedAmino = mutatedAminoMatcher.group(1).toUpperCase(Locale.ENGLISH);

          // Ignore terminating mutations, because they dont cause an aminoacid sequence alternation
          if ("TER".equals(mutatedAmino)) {
            continue;
          }

          String protein = sb.toString();

          boolean isDeletion = "DEL".equals(mutatedAmino);
          if (!isDeletion) {
            checkAminoAtMutationIndexMatches(protein, mutationIndex, mutatedAmino,
                header + " in line " + (headerLine + 1));
          }

          String output =
              getMutationSequence(protein, mutationIndex, params.getOffset(), isDeletion);
          String uniqueSolution = nativeAmino + ">" + mutatedAmino + ":" + output;
          if (uniqueSolutions.add(uniqueSolution)) {
            result.add(header);
            result.add(output);
          }
        }
      }
    }
    return result;
  }

  private static void checkAminoAtMutationIndexMatches(String protein, int mutationIndex,
      String expectedAmino, String mutationDescription) {
    char expectedAminoCharacter;
    if (expectedAmino.length() == 1) {
      expectedAminoCharacter = expectedAmino.charAt(0);
    } else {
      if (!AMINO_ACIDS.containsKey(expectedAmino)) {
        throw new IllegalArgumentException(
            String.format("Unknown aminoacid %s in mutation ", expectedAmino, mutationDescription));
      }
      expectedAminoCharacter = AMINO_ACIDS.get(expectedAmino);
    }

    char actualAminoCharacter = protein.charAt(mutationIndex - 1);
    if (expectedAminoCharacter != actualAminoCharacter) {
      throw new IllegalArgumentException(
          String.format("Incorrect aminoacid at index %s! Expected %s but got %s for mutation ",
              mutationIndex, expectedAminoCharacter, actualAminoCharacter, mutationDescription));
    }
  }

  /**
   * Retrieves the characters around mutationIndex (1 based) of the specified protein.
   * {@code offset} is the number of characters before and after the mutation that are retrieved.
   *
   * @param protein
   * @param mutationIndex (1 based)
   * @param offset the number of characters before and after the mutation that are retrieved
   * @return
   */
  private static String getMutationSequence(String protein, int mutationIndex, int offset,
      boolean isDeletion) {
    int zeroBasedIndex = mutationIndex - 1;
    int beginIndex = Math.max(0, zeroBasedIndex - offset);
    int endIndex = Math.min(protein.length(), zeroBasedIndex + (isDeletion ? 0 : 1) + offset);
    return protein.substring(beginIndex, endIndex);
  }
}
