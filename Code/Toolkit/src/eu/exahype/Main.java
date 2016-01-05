package eu.exahype;


public class Main {

  public static void main(String[] args) {
    if (args.length !=1 ) {
	  System.err.println("Please provide input file as argument" );
	  return;
    };

    String inputFileName = args[0];
    System.out.print("read input file " + inputFileName);

    eu.exahype.parser.Parser parser   = null;
    eu.exahype.node.Start    document = null;
	try {
		parser = new eu.exahype.parser.Parser(
		    new eu.exahype.lexer.Lexer(new java.io.PushbackReader(
		        new java.io.FileReader(inputFileName))));
	    document = parser.parse();
	    System.out.println(" ... ok");
	} catch (Exception e) {
	    System.out.println(" ... error: " + e.toString());
	}
  }
}
