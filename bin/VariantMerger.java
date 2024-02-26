//usr/bin/env jbang "$0" "$@" ; exit $?
//REPOS jcenter,genepi-maven=https://genepi.i-med.ac.at/maven
//DEPS info.picocli:picocli:4.6.1
//DEPS genepi:genepi-io:1.1.1

import java.io.File;
import java.util.concurrent.Callable;
import genepi.io.table.reader.CsvTableReader;
import genepi.io.table.writer.CsvTableWriter;
import picocli.CommandLine;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

public class VariantMerger implements Callable<Integer> {

	@Parameters(description = "Combined variants file")
	private String file;

	@Option(names = "--output", description = "Output files", required = true)
	private String output;

	@Option(names = "--indel-tag", description = "Detect indels by this tag", required = false)
	private String tag = "INDEL";

	public Integer call() throws Exception {

		assert (file != null);
		assert (output != null);

		CsvTableWriter writer = new CsvTableWriter(new File(output).getAbsolutePath(), '\t', false);
		CsvTableReader reader = new CsvTableReader(file, '\t');

		int lastPos = 0;
		String[] lastRow = null;
		int diff = 0;
		writer.setColumns(reader.getColumns());

		while (reader.next()) {

			int pos = reader.getInteger("Pos");
			int refLength = reader.getString("Ref").length();
			int variantlength = reader.getString("Variant").length();

			// init new position
			if (lastPos == 0) {
				diff = refLength - variantlength;
				lastPos = pos;
				lastRow = reader.getRow();
				continue;
			}

			if (comparePositions(lastPos, pos, diff)) {
				writer.setRow(lastRow);
				writer.next();
				lastPos = 0;
			} else {
				// since no hit, write previous row.
				writer.setRow(lastRow);
				writer.next();
				diff = refLength - variantlength;
				lastPos = pos;
				lastRow = reader.getRow();
			}

		}
		writer.setRow(lastRow);
		writer.next();
		reader.close();
		writer.close();
		return 0;
	}

	public static void main(String... args) {
		int exitCode = new CommandLine(new VariantMerger()).execute(args);
		System.exit(exitCode);
	}

	public boolean comparePositions(int lastPos, int pos, int diff) {
		for (int i = 0; i <= Math.abs(diff); i++) {
			if ((lastPos + i) == pos) {
				return true;
			}
		}
		return false;
	}
}
