package eu.exahype;


public class Main {

  public static void printHeader() {
	System.out.println( "================================" );
	System.out.println( " ___          _  _      ___ ___" ); 
	System.out.println( "/ __|_ ____ _| || |_  _/ _ \\ __|" );
	System.out.println( "| _|\\ \\ / _` | __ | || |  _/ _| " );
	System.out.println( "\\___/_\\_\\__,_|_||_|\\_, |_| \\___|" );
	System.out.println( "                   |__/         " );
	System.out.println( "================================" );
	System.out.println( "");
	System.out.println( " www.exahype.eu ");
	System.out.println( "");
	System.out.println( "================================" );	 
	System.out.println( "");
	System.out.println( "The project has received funding from the European Union's ");
	System.out.println( "Horizon 2020 research and innovation programme under grant ");
	System.out.println( "agreement No 671698 (ExaHyPE). It is based upon the PDE ");
	System.out.println( "framework Peano (www.peano-framework.org).");  
	System.out.println( "");
	System.out.println( "");
  }
  
  
  public static void waitForInteraction(boolean interactive) {
    if (interactive) {
	  System.out.println("<press Enter>");
	  try {
        System.in.read();
	  }
	  catch (Exception e) {}
	  for (int i = 0; i < 50; ++i) System.out.println();
	}
  }

  
  public static void main(String[] args) {
    printHeader();
    
    if (args.length !=1 && args.length !=2) {
	  System.err.println("ERROR: Please provide input file as argument" );
	  return;
    };
    
    if (args.length ==2 && args[0].compareTo("--not-interactive")!=0 && args[0].compareTo("--interactive")!=0 ) {
	  System.err.println("ERROR: First optional argument has to be --not-interactive or --interactive. Received \"" + args[0] + "\"" );
	  return;
    };
    
    boolean interactive = args.length==1 || args[0].compareTo("--interactive")==0;

    if (args.length==1) {
  	  System.out.println("INFO: You might want to add --not-interactive or --interactive as first command ");
  	  System.out.println("      line argument to control whether script runs interactively" );
    }
    	
    String inputFileName = args.length==2 ? args[1] : args[0];
    System.out.print("read input file " + inputFileName);

    eu.exahype.parser.Parser parser   = null;
    eu.exahype.node.Start    document = null;
	try {
	  parser = new eu.exahype.parser.Parser(
        new eu.exahype.lexer.Lexer(new java.io.PushbackReader(
		  new java.io.FileReader(inputFileName))));
	  document = parser.parse();
	  System.out.println(" ... ok");
	  System.out.println("Start to interpret script ... ");
	  waitForInteraction(interactive);
	} catch (Exception e) {
      System.out.println(" ... error: " + e.toString());
	  return;
	}
  }
}
