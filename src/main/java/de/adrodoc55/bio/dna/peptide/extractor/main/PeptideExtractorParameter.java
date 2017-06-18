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
import java.util.List;
import java.util.regex.Pattern;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;

/**
 * @author Adrodoc55
 */
public class PeptideExtractorParameter {
  @Parameter(names = {"-h", "--help"}, help = true,
      description = "Print information about the commandline usage")
  private boolean help;

  @Parameter(required = true, description = "<input-file>")
  private List<File> input;

  @Parameter(names = {"-o", "--output"}, required = true, description = "Specify an output file")
  private File output;

  @Parameter(names = {"--offset"},
      description = "The number of aminoacids before and after each mutation that are retrieved")
  private int offset = 8;

  @Parameter(names = {"-p", "--header-prefix"},
      description = "The prefix to identify a mutation header line")
  private String headerPrefix = ">";

  @Parameter(names = {"-r", "--header-regex"},
      description = "The regex to match mutation header lines (header lines that don't match will be ignored)")
  private String headerRegex = ">.*[A-Za-z]{3}\\d+[A-Za-z]{3}$";

  @Parameter(names = {"-i", "--mutation-index-regex"},
      description = "The regex for retrieving the mutation index (1 based) from a header line (group 1 will be used)")
  private String mutationIndexRegex = ">.*[A-Za-z]{3}(\\d+)[A-Za-z]{3}$";

  @Parameter(names = {"-n", "--native-amino-regex"},
      description = "The regex for retrieving the native aminoacid from a header line (group 1 will be used)")
  private String nativeAminoRegex = ">.*([A-Za-z]{3})\\d+[A-Za-z]{3}$";

  @Parameter(names = {"-m", "--mutated-amino-regex"},
      description = "The regex for retrieving the mutated aminoacid from a header line (group 1 will be used)")
  private String mutatedAminoRegex = ">.*[A-Za-z]{3}\\d+([A-Za-z]{3})$";

  public boolean isHelp() {
    return help;
  }

  public File getInput() throws ParameterException {
    if (input.size() != 1) {
      throw new ParameterException("Exactly one source file has to be specified");
    }
    return input.get(0).getAbsoluteFile();
  }

  public File getOutput() {
    return output;
  }

  public int getOffset() {
    return offset;
  }

  public String getHeaderPrefix() {
    return headerPrefix;
  }

  private Pattern headerPattern;

  public Pattern getHeaderPattern() {
    if (headerPattern == null) {
      headerPattern = Pattern.compile(headerRegex);
    }
    return headerPattern;
  }

  private Pattern mutationIndexPattern;

  public Pattern getMutationIndexPattern() {
    if (mutationIndexPattern == null) {
      mutationIndexPattern = Pattern.compile(mutationIndexRegex);
    }
    return mutationIndexPattern;
  }

  private Pattern nativeAminoPattern;

  public Pattern getNativeAminoPattern() {
    if (nativeAminoPattern == null) {
      nativeAminoPattern = Pattern.compile(nativeAminoRegex);
    }
    return nativeAminoPattern;
  }

  private Pattern mutatedAminoPattern;

  public Pattern getMutatedAminoPattern() {
    if (mutatedAminoPattern == null) {
      mutatedAminoPattern = Pattern.compile(mutatedAminoRegex);
    }
    return mutatedAminoPattern;
  }
}
