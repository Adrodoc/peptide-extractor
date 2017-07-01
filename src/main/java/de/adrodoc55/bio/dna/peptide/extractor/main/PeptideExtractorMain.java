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

import static de.adrodoc55.bio.dna.FastaConstants.COMMENT_PREFIX;
import static de.adrodoc55.bio.dna.FastaConstants.HEADER_PREFIX;
import static de.adrodoc55.bio.dna.FastaConstants.MAX_LINE_LENGTH;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.HashSet;
import java.util.Set;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParameterException;
import com.google.common.base.Splitter;
import com.google.common.io.Files;

import de.adrodoc55.bio.dna.peptide.extractor.PeptideExtractorException;
import de.adrodoc55.bio.dna.peptide.extractor.mutation.Mutation;
import de.adrodoc55.bio.dna.peptide.extractor.mutation.Mutations;
import de.adrodoc55.bio.dna.peptide.extractor.mutation.SilentSingleNucleotidePolymorphism;
import de.adrodoc55.bio.dna.peptide.extractor.mutation.Termination;

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

      main(params);
    } catch (ParameterException ex) {
      System.err.println(ex.getLocalizedMessage());
      System.err.println("Run with '-h' to print help");
    }
  }

  private static void main(PeptideExtractorParameter params)
      throws IOException, PeptideExtractorException {
    File input = params.getInput();

    File outputFile = params.getOutput();
    Files.createParentDirs(outputFile);
    outputFile.createNewFile();

    try (LineNumberReader in = new LineNumberReader(new FileReader(input));
        BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));) {
      Set<String> uniqueSolutions = new HashSet<>();
      String lastHeader = null;
      int lastHeaderLineNumber = -1;
      StringBuilder protein = new StringBuilder();
      while (true) {
        String line = in.readLine();
        if (line == null || line.startsWith(HEADER_PREFIX)) {
          if (lastHeader != null) {
            processMutation(lastHeader, lastHeaderLineNumber, protein, uniqueSolutions, out,
                params);
          }
        } else if (!line.startsWith(COMMENT_PREFIX)) {
          protein.append(line.trim());
        }
        if (line == null) {
          break;
        } else if (line.startsWith(HEADER_PREFIX)) {
          lastHeader = line;
          lastHeaderLineNumber = in.getLineNumber();
          protein.setLength(0);
        }
      }
    }
  }

  private static final Splitter SPLITTER = Splitter.fixedLength(MAX_LINE_LENGTH);

  private static void processMutation(String header, int headerLineNumber, StringBuilder protein,
      Set<String> uniqueSolutions, BufferedWriter out, PeptideExtractorParameter params)
      throws IOException, PeptideExtractorException {
    String mutationDescription = header + " in line " + headerLineNumber;
    try {
      Mutation mutation = Mutations.parse(header);
      if (mutation instanceof SilentSingleNucleotidePolymorphism) {
        System.err.println("Ignoring silent mutation " + mutationDescription
            + ", because silent mutations don't cause an aminoacid sequence alternation");
      } else if (mutation instanceof Termination) {
        System.err.println("Ignoring terminating mutation " + mutationDescription
            + ", because terminating mutations don't cause an aminoacid sequence alternation");
      } else if (mutation != null) {
        CharSequence output = mutation.extractFromProtein(protein, params.getEnclosing());
        String uniqueSolution = mutation.getUniqueSolution(output);
        if (uniqueSolutions.add(uniqueSolution)) {
          out.write(header);
          out.newLine();
          for (String subOutput : SPLITTER.split(output)) {
            out.write(subOutput);
            out.newLine();
          }
        }
      } else {
        throw new PeptideExtractorException("Unrecognized mutation header");
      }
    } catch (PeptideExtractorException ex) {
      if (params.isIgnoreErrors()) {
        System.err.println(
            "Ignoring mutation " + mutationDescription + " due to: " + ex.getLocalizedMessage());
      } else {
        throw new PeptideExtractorException(
            "Error at mutation " + mutationDescription + ": " + ex.getLocalizedMessage(), ex);
      }
    }
  }
}
