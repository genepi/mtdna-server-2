//usr/bin/env jbang "$0" "$@" ; exit $?
//REPOS jcenter,genepi-maven=https://genepi.i-med.ac.at/maven
//DEPS info.picocli:picocli:4.6.1
//DEPS genepi:genepi-io:1.1.1

import java.io.File;
import java.util.concurrent.Callable;

import org.apache.commons.math3.distribution.BinomialDistribution;

import genepi.io.table.reader.CsvTableReader;
import genepi.io.table.writer.CsvTableWriter;
import picocli.CommandLine;
import picocli.CommandLine.Option;

public class CoverageCorrection implements Callable<Integer> {

	@Option(names = "--input", description = "Input Variants File", required = true)
	private String input;

	@Option(names = "--output", description = "Output Variants File", required = true)
	private String output;

	public Integer call() throws Exception {

		CsvTableReader reader = new CsvTableReader(input, '\t');
		CsvTableWriter writer = new CsvTableWriter(new File(output).getAbsolutePath(), '\t', false);
		writer.setColumns(reader.getColumns());

		while (reader.next()) {
			String[] currentRow = reader.getRow();
			if (reader.getString("Coverage").contains(".") || reader.getString("Coverage").equals("0")) {
				writer.setRow(reader.getRow());
				writer.next();
				continue;
			}
			if (reader.getString("VariantLevel").contains(",")) {
				writer.setRow(reader.getRow());
				writer.next();
				continue;
			}
			
			if (reader.getString("MeanBaseQuality").contains(".")) {
				writer.setRow(reader.getRow());
				writer.next();
				continue;
			}
			
			int coverage = reader.getInteger("Coverage");
			double variantLevel = reader.getDouble("VariantLevel");
			double meanBaseQuality = reader.getDouble("MeanBaseQuality");

			double error = phredToProbability(meanBaseQuality);
			double errLim = 0.999;
			double detLim = 0.999;

			double estimatedLevel = relaxTPFPMin(coverage, error, errLim, detLim, 1);

			if (estimatedLevel > variantLevel) {
				currentRow[reader.getColumnIndex("Filter")] = "LOW_COVERAGE";
			}
			writer.setRow(currentRow);

			writer.next();
		}
		writer.close();
		return 0;
	}

	public static void main(String... args) {
		int exitCode = new CommandLine(new CoverageCorrection()).execute(args);
		System.exit(exitCode);
	}

	public static final int INITIAL_COVERAGE = 1;

	public static double relaxTPFPMin(int coverage, double err, double errLimit, double detLimit, int minVar) {
		double resVaf = 0;
		double currentVaf = 0.001; // Start with a very low VAF

		while (true) {
			// Relaxation condition no 1.
			int errs = 0;
			while (new BinomialDistribution(coverage, err).cumulativeProbability(errs) < errLimit) {
				errs++;
			}

			// Relaxation condition no 2.
			int variants = errs;
			if (errs == 0) {
				variants = 1;
			}

			double aux = 1 - new BinomialDistribution(coverage, currentVaf).cumulativeProbability(variants - 1);

			if (aux >= detLimit && minVar <= variants) {
				resVaf = currentVaf;
				break;
			} else {
				currentVaf += 0.001; // Increment VAF by a small step
			}
		}

		return resVaf; // Return minimum VAF and coverage
	}

	public static double phredToProbability(double meanBaseQuality) {
		return Math.pow(10, -((double) meanBaseQuality / 10));
	}

}
