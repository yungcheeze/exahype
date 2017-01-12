package eu.exahype;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Iterator;

/**
 * This method implements a mechanism to allow files to be at any subdirectory
 * from they are expected. This allows advanced code structuring. For instance:
 *    File a = new File("basedir/foo.c");
 * In order to structure code, this file was moved after toolkit running to
 *    basedir/FooStuff/foo.c
 * It is important that the filename does not change and is sufficiently unique
 * to allow the following:
 *    a = FileSearch.relocate(a);
 * This will return a new File("basedir/FooStuff/foo.c") and print out information
 * at System.out.println.
 *
 * If you find this to be a performance bottleneck for large directory trees,
 * consider rewriting it with a <b>breadth-first search traversal method which quits
 * after the first occurence of a result</b>.
 * This proposed algorithm implements the idea that an ExaHyPE-related file is
 * typically not nested deeply in an (eg.) large external library in the same
 * directory tree.
 *
 **/
public class FileSearch {
	/**
	 * Returns true if there is another candidate for file `original`.
	 **/
	public static boolean relocatable(File original) {
		File replacement = relocate(original);
		return !(original.getAbsolutePath() == replacement.getAbsolutePath());
	}

	/**
	 * Tries to relocate a file. If not relocatable, it just returns the original
	 * file. In the spirit of the ExaHyPE toolkit, it just writes out to stdout if
	 * it has the feeling to do so.
	 **/
	public static File relocate(File original) {
		File basedir = original.getParentFile();
		String filename = original.getName();
		
		// basedir does not exist, so no place to look in.
		if(!basedir.isDirectory())
			return original;
		
		List<String> results = new FileSearch(basedir.getAbsolutePath(), filename).getResult();
		
		if(results.size()==0)
			return original;

		// filter out equal items
		Iterator<String> it = results.iterator();
		while(it.hasNext()) {
			String candidate = it.next();
			if( candidate == original.getAbsolutePath() )
				it.remove();
		}
		
		if(results.size()==0)
			return original;

		if(results.size()==1) {
			System.out.println("Found candidate for "+filename+" at "+results.get(0));
			return new File(results.get(0));
		}

		if(results.size()>1) {
			System.out.println("Found multiple candidates for '" + original.getAbsolutePath() + "':");
			for(int i = 0; i < results.size(); i++)
				System.out.println(results.get(i));
			System.out.println("Sticking to original file.");
			return original;
		}
		
		return original;
	}
	
	/**
	 * Syntactic sugar which allows you to replace `new java.io.File(...)` expressions
	 * by `FileSearch.relocatableFile(...)` without the need of an additional line
	 * as in
	 *
	 *   File a = new java.io.File(...);
	 *   a = FileSearch.relocate(a);
	 *
	 **/
	public static File relocatableFile(String original) {
		return relocate(new File(original));
	}

	/* Java resursive file finding. Ugly as hell. */
	
	private String fileNameToSearch;
	private List<String> result = new ArrayList<String>();
	public List<String> getResult() { return result; }
	
	public FileSearch(String basedir, String filename) {
		fileNameToSearch = filename;
		search(new File(basedir));
	}

	public void search(File file) {
		if (file.isDirectory()) {
			if (file.canRead()) {
				for (File temp : file.listFiles()) {
					if (temp.isDirectory()) {
						search(temp);
					} else {
						if(fileNameToSearch.equals(temp.getName()))
							result.add(temp.getAbsolutePath());
					}
				}
			}
		}
	}
	

	/**
	 * For ad hoc unit testing, use the class as standalone CLI program.
	 *
	 * Compile like:
	 *   javac eu/exahype/FileSearch.java
	 *
	 * Run like:
	 *    mkdir -p /tmp/testroot/a/b/c
	 *    touch /tmp/testroot/a/b/c/Demo.txt
	 *    java eu.exahype.FileSearch  /tmp/testroot/Demo.txt
	 *
	 * Output is:
	 *
	 *    Start: /tmp/testroot/Demo.txt
	 *    End: /tmp/testroot/a/b/c/Demo.txt
	 *
	 * If you add another file, it sticks to the original:
	 *
	 *     touch /tmp/testroot/a/b/Demo.txt
	 *     java eu.exahype.FileSearch  /tmp/testroot/Demo.txt
	 *
	 * Output then is:
	 *
	 *     Found multiple candidates for '/tmp/testroot/Demo.txt':
	 *     /tmp/testroot/a/b/c/Demo.txt
	 *     /tmp/testroot/a/b/Demo.txt
	 *     Sticking to original file.
	 *     Start: /tmp/testroot/Demo.txt
	 *     End: /tmp/testroot/Demo.txt
	 *
	 * So it apparently works.
	 */
	public static void main(String[] args) {
		if(args.length != 1 || args[0] == "--help") {
			System.out.println("Usage: java FileSearch your/file/name.c");
			System.exit(1);
		}
		
		File start = new File(args[0]);
		File end = relocate(start);
		
		System.out.println("Start: "+start.getAbsolutePath());
		System.out.println("End: "+end.getAbsolutePath());
	}
}