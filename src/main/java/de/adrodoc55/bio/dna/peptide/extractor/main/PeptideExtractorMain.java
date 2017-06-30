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

import static de.adrodoc55.bio.dna.FastaConstants.HEADER_PREFIX;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParameterException;
import com.google.common.base.Charsets;
import com.google.common.base.Joiner;
import com.google.common.collect.ImmutableList;
import com.google.common.io.Files;
import com.google.common.io.Resources;

import de.adrodoc55.bio.dna.peptide.extractor.PeptideExtractorException;
import de.adrodoc55.bio.dna.peptide.extractor.mutation.Mutation;
import de.adrodoc55.bio.dna.peptide.extractor.mutation.Mutations;

/**
 * @author Adrodoc55
 */
public class PeptideExtractorMain {
  public static void main(String[] args) throws IOException, PeptideExtractorException {
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

  private static List<String> transform(List<String> lines, PeptideExtractorParameter params)
      throws PeptideExtractorException {
    List<String> result = new ArrayList<>();
    Set<String> uniqueSolutions = new HashSet<>();
    for (int lineIndex = 0; lineIndex < lines.size();) {
      int headerLineIndex = lineIndex++;
      String header = lines.get(headerLineIndex);
      if (header.startsWith(HEADER_PREFIX)) {

        StringBuilder sb = new StringBuilder();
        while (lineIndex < lines.size() && !lines.get(lineIndex).startsWith(HEADER_PREFIX)) {
          sb.append(lines.get(lineIndex++));
        }

        int headerLineNumber = headerLineIndex + 1;
        try {
          Mutation mutation = Mutations.parse(header);
          if (mutation != null) {
            String protein = sb.toString();
            String output = mutation.extractFromProtein(protein, params.getOffset());
            String uniqueSolution = mutation.getUniqueSolution(output);
            if (uniqueSolutions.add(uniqueSolution)) {
              result.add(header);
              result.add(output);
            }
          } else {
            throw new PeptideExtractorException("Unrecognized mutation header");
          }
        } catch (PeptideExtractorException ex) {
          String mutationDescription = header + " in line " + headerLineNumber;
          if (params.isIgnoreErrors()) {
            System.err.println("Ignoring mutation " + mutationDescription + " due to: "
                + ex.getLocalizedMessage());
          } else {
            throw new PeptideExtractorException(
                "Error at mutation " + mutationDescription + ": " + ex.getLocalizedMessage(), ex);
          }
        }
      }
    }
    return result;
  }
}
