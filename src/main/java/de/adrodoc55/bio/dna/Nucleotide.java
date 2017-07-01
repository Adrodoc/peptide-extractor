package de.adrodoc55.bio.dna;

import java.util.Locale;

import de.adrodoc55.bio.dna.peptide.extractor.UnknownNucleotideException;

public enum Nucleotide {
  A, C, G, T;
  public static Nucleotide fromCharCode(char code) throws UnknownNucleotideException {
    String upperCase = String.valueOf(code).toUpperCase(Locale.ENGLISH);
    try {
      return valueOf(upperCase);
    } catch (IllegalArgumentException ex) {
      throw new UnknownNucleotideException("Unknown nucleotide " + code, ex);
    }
  }
}
