package meshi.applications.TriC.yeast;

import java.io.BufferedWriter;
import java.io.FileWriter;

import meshi.applications.TriC.TricYeastAlignment;
import meshi.applications.prediction.GDTcalculator;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Residue;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.MeshiProgram;
import meshi.util.file.File2StringArray;

public class Willison64PossibilitiesYaoLikeSingle extends MeshiProgram implements Residues,AtomTypes {
	
	public static void main(String[] args) throws Exception {
		init(args);
		
		// Defining the arrangement 
		String topTrue = "ZQHEBDAG"; // OMS
//		String topTrue = "ZQHEBDAG"; // OMS+1
//		String topTrue = "ZQHEBDAG"; // OMS-1
//		String topTrue = "ZEAHDQGB"; // Willison
//		String topTrue = "AZBGQDEH"; // Yao
		String chainsTop = "FEAGDHCB";
		String chainsTop_E = "feagdhcb";
		String genesTop =  "ZEAHDQGB";

		String botTrue = "ZQHEBDAG"; // OMS
//		String botTrue = "GZQHEBDA"; // OMS+1
//		String botTrue = "QHEBDAGZ"; // OMS-1
//		String botTrue = "ZEAHDQGB"; // Willison
//		String botTrue = "AZBGQDEH"; // Yao
		String chainsBot = "NMIOLPKJ";
		String chainsBot_E = "nmiolpkj";
		String genesBot =  "ZEAHDQGB";

		String botToWrite = ""+botTrue.charAt(0);
		for (int c=7 ; c>0 ; c--) {
			botToWrite = botToWrite + botTrue.charAt(c);
		}

		AtomList willison1Q3R_D = turnWillisonTo1Q3R(new AtomList("/Users/nirka/TRiC/New_Work_on_TRiC_14_8_2011/Gunnar/3P9D.pdb"),
				chainsTop, genesTop, chainsBot, genesBot,false);
		Atom.resetNumberOfAtoms();
		AtomList willison1Q3R_E = turnWillisonTo1Q3R(new AtomList("/Users/nirka/TRiC/New_Work_on_TRiC_14_8_2011/Gunnar/3P9E.pdb"),
				chainsTop_E, genesTop, chainsBot_E, genesBot,true);
		Atom.resetNumberOfAtoms();
		AtomList unit = new AtomList("/Users/nirka/TRiC/New_Work_on_TRiC_14_8_2011/Gunnar/1Q3R_A.pdb");					
		AtomList equaUnit = filterDomainAccordingTo1Q3R(unit, 'K' , 'E');
		AtomList middleUnit = filterDomainAccordingTo1Q3R(unit, 'K' , 'M');
		AtomList apicUnit = filterDomainAccordingTo1Q3R(unit, 'K' , 'A');
		TricYeastAlignment alignments = new TricYeastAlignment();

		// The chain2gene files:
		int modelNumber = 0;
		for (int delta1=0 ; delta1<8 ; delta1++) { // moving the real 2-fold "delta1" units in the top ring order from willisons 2-fold. In complex 3P9E
			for (int delta2=0 ; delta2<8 ; delta2++) { // moving the real 2-fold "delta1" units in the top ring order from willisons 2-fold. In complex 3P9E
				String outputFile;
				if (modelNumber>9) {
					outputFile = "chain2gene_"+modelNumber+".txt";
				}
				else {
					outputFile = "chain2gene_0"+modelNumber+".txt";
				}
				BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile));
				modelNumber++;
				// Doing Willison 3P9D
				// Doing top ring
				for (int pos1=0 ; pos1<8 ; pos1++) { // Position on Willison's PDB
					char chainID = chainsTop.charAt(pos1);
					char toPut = topTrue.charAt((pos1+8-delta1) % 8);
					bw.write(chainID + " " + toPut + "\n");
				}
				// Doing bottom ring
				for (int pos1=0 ; pos1<8 ; pos1++) { // Position on Willison's PDB
					char chainID = chainsBot.charAt(pos1);
					char toPut = botTrue.charAt((pos1+delta1) % 8);
					bw.write(chainID + " " + toPut + "\n");
				}
				// Doing Willison 3P9E
				// Doing top ring
				for (int pos1=0 ; pos1<8 ; pos1++) { // Position on Willison's PDB
					char chainID = chainsTop_E.charAt(pos1);
					char toPut = topTrue.charAt((pos1+8-delta2) % 8);
					bw.write(chainID + " " + toPut + "\n");
				}
				// Doing bottom ring
				for (int pos1=0 ; pos1<8 ; pos1++) { // Position on Willison's PDB
					char chainID = chainsBot_E.charAt(pos1);
					char toPut = botTrue.charAt((pos1+delta2) % 8);
					bw.write(chainID + " " + toPut + "\n");
				}
				bw.close();
			}		
		}

		
//		try {
		modelNumber = 0;
		for (int delta1=0 ; delta1<8 ; delta1++) { // moving the real 2-fold "delta1" units in the top ring order from willisons 2-fold. In complex 3P9E
			for (int delta2=0 ; delta2<8 ; delta2++) { // moving the real 2-fold "delta1" units in the top ring order from willisons 2-fold. In complex 3P9E
					String outputFile;
					if (modelNumber>9) {
						outputFile = "model_"+modelNumber+"_"+topTrue+"_"+botTrue+"_"+delta1+"_"+delta2+".pdb";
					}
					else {
						outputFile = "model_0"+modelNumber+"_"+topTrue+"_"+botTrue+"_"+delta1+"_"+delta2+".pdb";
					}
					BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile));
					modelNumber++;
					// Doing Willison 3P9D
					// Doing top ring
					for (int pos1=0 ; pos1<8 ; pos1++) { // Position on Willison's PDB
						char chainID = chainsTop.charAt(pos1);
						char toPut = topTrue.charAt((pos1+8-delta1) % 8);
						AtomList template = willison1Q3R_D.chainFilter(chainID+"");
						// Doing Equatiroial
						GDTcalculator.alignBySubset(template, equaUnit, 0.75);
						// Doing Middle
						GDTcalculator.alignBySubset(template, middleUnit, 0.75);
						// Doing Apical
						GDTcalculator.alignBySubset(template, apicUnit, 0.75);
						// Finally doing the alignment 
						Atom.resetNumberOfAtoms();
						AtomList finalUnit = HMunit(alignments.getAlignment(toPut+""), toPut,alignments.getAlignment("K"), 'K',unit);
						finalUnit.setChain(chainID+"");
						// Writing
						for (int c=0 ; c<finalUnit.size() ; c++) {
							bw.write(finalUnit.atomAt(c).toString() + "\n");
						}
						bw.write("TER\n");			
					}
					// Doing bottom ring
					for (int pos1=0 ; pos1<8 ; pos1++) { // Position on Willison's PDB
						char chainID = chainsBot.charAt(pos1);
						char toPut = botTrue.charAt((pos1+delta1) % 8);
						AtomList template = willison1Q3R_D.chainFilter(chainID+"");
						// Doing Equatiroial
						GDTcalculator.alignBySubset(template, equaUnit, 0.75);
						// Doing Middle
						GDTcalculator.alignBySubset(template, middleUnit, 0.75);
						// Doing Apical
						GDTcalculator.alignBySubset(template, apicUnit, 0.75);
						// Finally doing the alignment 
						Atom.resetNumberOfAtoms();
						AtomList finalUnit = HMunit(alignments.getAlignment(toPut+""), toPut,alignments.getAlignment("K"), 'K',unit);
						finalUnit.setChain(chainID+"");
						// Writing
						for (int c=0 ; c<finalUnit.size() ; c++) {
							bw.write(finalUnit.atomAt(c).toString() + "\n");
						}
						bw.write("TER\n");			
					}
					// Doing Willison 3P9E
					// Doing top ring
					for (int pos1=0 ; pos1<8 ; pos1++) { // Position on Willison's PDB
						char chainID = chainsTop_E.charAt(pos1);
						char toPut = topTrue.charAt((pos1+8-delta2) % 8);
						AtomList template = willison1Q3R_E.chainFilter(chainID+"");
						// Doing Equatiroial
						GDTcalculator.alignBySubset(template, equaUnit, 0.75);
						// Doing Middle
						GDTcalculator.alignBySubset(template, middleUnit, 0.75);
						// Doing Apical
						GDTcalculator.alignBySubset(template, apicUnit, 0.75);
						// Finally doing the alignment 
						Atom.resetNumberOfAtoms();
						AtomList finalUnit = HMunit(alignments.getAlignment(toPut+""), toPut,alignments.getAlignment("K"), 'K',unit);
						finalUnit.setChain(chainID+"");
						// Writing
						for (int c=0 ; c<finalUnit.size() ; c++) {
							bw.write(finalUnit.atomAt(c).toString() + "\n");
						}
						bw.write("TER\n");			
					}
					// Doing bottom ring
					for (int pos1=0 ; pos1<8 ; pos1++) { // Position on Willison's PDB
						char chainID = chainsBot_E.charAt(pos1);
						char toPut = botTrue.charAt((pos1+delta2) % 8);
						AtomList template = willison1Q3R_E.chainFilter(chainID+"");
						// Doing Equatiroial
						GDTcalculator.alignBySubset(template, equaUnit, 0.75);
						// Doing Middle
						GDTcalculator.alignBySubset(template, middleUnit, 0.75);
						// Doing Apical
						GDTcalculator.alignBySubset(template, apicUnit, 0.75);
						// Finally doing the alignment 
						Atom.resetNumberOfAtoms();
						AtomList finalUnit = HMunit(alignments.getAlignment(toPut+""), toPut,  alignments.getAlignment("K"), 'K', unit);
						finalUnit.setChain(chainID+"");
						// Writing
						for (int c=0 ; c<finalUnit.size() ; c++) {
							bw.write(finalUnit.atomAt(c).toString() + "\n");
						}
						bw.write("TER\n");			
					}
					bw.write("END\n");			
					bw.close();
				}		
			}
//		}
//		catch(Exception e) {
//			throw new RuntimeException(e.getMessage());
//		}

		
	}

	
	protected static AtomList HMunit(String querySeq,char queryUnit , String templateSeq, char templateUnit, AtomList template) {
		boolean[] residueInAllSeqs = {false,  //  0
				false,  // 1 
				false,  // 2 
				false,  // 3 
				false,  // 4 
				false,  // 5 
				false,  // 6 
				false,  // 7 
				false,  // 8 
				false,  // 9 
				false,  // 10 
				false,  // 11 
				false,  // 12 
				false,  // 13 
				false,  // 14 
				false,  // 15 
				false,  // 16 
				false,  // 17 
				false,  // 18 
				false,  // 19 
				false,  // 20 
				false,  // 21 
				false,  // 22 
				false,  // 23 
				false,  // 24 
				false,  // 25 
				false,  // 26 
				false,  // 27 
				true,  // 28 
				true,  // 29 
				true,  // 30 
				true,  // 31 
				true,  // 32 
				true,  // 33 
				true,  // 34 
				true,  // 35 
				true,  // 36 
				true,  // 37 
				true,  // 38 
				true,  // 39 
				true,  // 40 
				true,  // 41 
				true,  // 42 
				true,  // 43 
				true,  // 44 
				true,  // 45 
				true,  // 46 
				true,  // 47 
				true,  // 48 
				true,  // 49 
				true,  // 50 
				true,  // 51 
				true,  // 52 
				true,  // 53 
				true,  // 54 
				true,  // 55 
				true,  // 56 
				true,  // 57 
				true,  // 58 
				true,  // 59 
				true,  // 60 
				true,  // 61 
				true,  // 62 
				true,  // 63 
				true,  // 64 
				true,  // 65 
				true,  // 66 
				true,  // 67 
				true,  // 68 
				true,  // 69 
				true,  // 70 
				true,  // 71 
				true,  // 72 
				true,  // 73 
				true,  // 74 
				true,  // 75 
				true,  // 76 
				false,  // 77 
				false,  // 78 
				true,  // 79 
				true,  // 80 
				true,  // 81 
				true,  // 82 
				true,  // 83 
				true,  // 84 
				true,  // 85 
				true,  // 86 
				true,  // 87 
				true,  // 88 
				true,  // 89 
				true,  // 90 
				true,  // 91 
				true,  // 92 
				true,  // 93 
				true,  // 94 
				true,  // 95 
				true,  // 96 
				true,  // 97 
				true,  // 98 
				true,  // 99 
				true,  // 100 
				true,  // 101 
				true,  // 102 
				true,  // 103 
				true,  // 104 
				true,  // 105 
				true,  // 106 
				true,  // 107 
				true,  // 108 
				true,  // 109 
				true,  // 110 
				true,  // 111 
				true,  // 112 
				true,  // 113 
				true,  // 114 
				true,  // 115 
				true,  // 116 
				true,  // 117 
				true,  // 118 
				true,  // 119 
				true,  // 120 
				true,  // 121 
				true,  // 122 
				true,  // 123 
				true,  // 124 
				true,  // 125 
				true,  // 126 
				true,  // 127 
				true,  // 128 
				true,  // 129 
				true,  // 130 
				true,  // 131 
				true,  // 132 
				true,  // 133 
				true,  // 134 
				true,  // 135 
				true,  // 136 
				true,  // 137 
				true,  // 138 
				false,  // 139 
				true,  // 140 
				true,  // 141 
				true,  // 142 
				true,  // 143 
				true,  // 144 
				true,  // 145 
				true,  // 146 
				true,  // 147 
				true,  // 148 
				true,  // 149 
				true,  // 150 
				true,  // 151 
				true,  // 152 
				true,  // 153 
				true,  // 154 
				true,  // 155 
				true,  // 156 
				true,  // 157 
				true,  // 158 
				true,  // 159 
				true,  // 160 
				true,  // 161 
				true,  // 162 
				true,  // 163 
				true,  // 164 
				true,  // 165 
				true,  // 166 
				true,  // 167 
				true,  // 168 
				true,  // 169 
				true,  // 170 
				true,  // 171 
				false,  // 172 
				false,  // 173 
				false,  // 174 
				false,  // 175 
				true,  // 176 
				true,  // 177 
				true,  // 178 
				true,  // 179 
				true,  // 180 
				true,  // 181 
				true,  // 182 
				true,  // 183 
				true,  // 184 
				true,  // 185 
				true,  // 186 
				true,  // 187 
				true,  // 188 
				true,  // 189 
				true,  // 190 
				true,  // 191 
				true,  // 192 
				true,  // 193 
				true,  // 194 
				true,  // 195 
				true,  // 196 
				true,  // 197 
				true,  // 198 
				true,  // 199 
				true,  // 200 
				true,  // 201 
				true,  // 202 
				true,  // 203 
				true,  // 204 
				true,  // 205 
				true,  // 206 
				true,  // 207 
				true,  // 208 
				true,  // 209 
				true,  // 210 
				true,  // 211 
				true,  // 212 
				true,  // 213 
				true,  // 214 
				true,  // 215 
				false,  // 216 
				false,  // 217 
				false,  // 218 
				false,  // 219 
				false,  // 220 
				false,  // 221 
				false,  // 222 
				false,  // 223 
				false,  // 224 
				false,  // 225 
				false,  // 226 
				false,  // 227 
				false,  // 228 
				true,  // 229 
				true,  // 230 
				true,  // 231 
				true,  // 232 
				true,  // 233 
				true,  // 234 
				true,  // 235 
				true,  // 236 
				true,  // 237 
				true,  // 238 
				true,  // 239 
				true,  // 240 
				true,  // 241 
				true,  // 242 
				true,  // 243 
				true,  // 244 
				true,  // 245 
				true,  // 246 
				true,  // 247 
				true,  // 248 
				true,  // 249 
				true,  // 250 
				true,  // 251 
				true,  // 252 
				true,  // 253 
				true,  // 254 
				true,  // 255 
				true,  // 256 
				true,  // 257 
				false,  // 258 
				false,  // 259 
				false,  // 260 
				false,  // 261 
				false,  // 262 
				true,  // 263 
				true,  // 264 
				true,  // 265 
				true,  // 266 
				true,  // 267 
				false,  // 268 
				false,  // 269 
				true,  // 270 
				true,  // 271 
				true,  // 272 
				true,  // 273 
				true,  // 274 
				true,  // 275 
				false,  // 276 
				false,  // 277 
				false,  // 278 
				false,  // 279 
				false,  // 280 
				false,  // 281 
				true,  // 282 
				true,  // 283 
				true,  // 284 
				true,  // 285 
				true,  // 286 
				true,  // 287 
				true,  // 288 
				true,  // 289 
				true,  // 290 
				true,  // 291 
				true,  // 292 
				true,  // 293 
				true,  // 294 
				true,  // 295 
				true,  // 296 
				true,  // 297 
				true,  // 298 
				false,  // 299 
				true,  // 300 
				true,  // 301 
				true,  // 302 
				true,  // 303 
				true,  // 304 
				true,  // 305 
				true,  // 306 
				true,  // 307 
				true,  // 308 
				true,  // 309 
				true,  // 310 
				true,  // 311 
				true,  // 312 
				true,  // 313 
				true,  // 314 
				true,  // 315 
				true,  // 316 
				true,  // 317 
				true,  // 318 
				true,  // 319 
				true,  // 320 
				true,  // 321 
				true,  // 322 
				true,  // 323 
				true,  // 324 
				true,  // 325 
				true,  // 326 
				true,  // 327 
				true,  // 328 
				true,  // 329 
				true,  // 330 
				true,  // 331 
				true,  // 332 
				true,  // 333 
				true,  // 334 
				true,  // 335 
				true,  // 336 
				true,  // 337 
				true,  // 338 
				true,  // 339 
				true,  // 340 
				true,  // 341 
				false,  // 342 
				false,  // 343 
				false,  // 344 
				false,  // 345 
				false,  // 346 
				false,  // 347 
				false,  // 348 
				false,  // 349 
				false,  // 350 
				false,  // 351 
				true,  // 352 
				true,  // 353 
				true,  // 354 
				true,  // 355 
				true,  // 356 
				true,  // 357 
				true,  // 358 
				true,  // 359 
				true,  // 360 
				true,  // 361 
				true,  // 362 
				true,  // 363 
				true,  // 364 
				true,  // 365 
				true,  // 366 
				true,  // 367 
				true,  // 368 
				true,  // 369 
				true,  // 370 
				true,  // 371 
				true,  // 372 
				true,  // 373 
				true,  // 374 
				true,  // 375 
				true,  // 376 
				true,  // 377 
				true,  // 378 
				true,  // 379 
				true,  // 380 
				true,  // 381 
				true,  // 382 
				true,  // 383 
				true,  // 384 
				true,  // 385 
				true,  // 386 
				true,  // 387 
				true,  // 388 
				true,  // 389 
				true,  // 390 
				true,  // 391 
				true,  // 392 
				true,  // 393 
				true,  // 394 
				true,  // 395 
				true,  // 396 
				true,  // 397 
				false,  // 398 
				false,  // 399 
				false,  // 400 
				false,  // 401 
				false,  // 402 
				false,  // 403 
				true,  // 404 
				true,  // 405 
				true,  // 406 
				true,  // 407 
				true,  // 408 
				true,  // 409 
				true,  // 410 
				true,  // 411 
				false,  // 412 
				true,  // 413 
				true,  // 414 
				true,  // 415 
				true,  // 416 
				true,  // 417 
				true,  // 418 
				true,  // 419 
				true,  // 420 
				false,  // 421 
				false,  // 422 
				true,  // 423 
				true,  // 424 
				true,  // 425 
				true,  // 426 
				true,  // 427 
				true,  // 428 
				true,  // 429 
				true,  // 430 
				true,  // 431 
				true,  // 432 
				true,  // 433 
				true,  // 434 
				false,  // 435 
				false,  // 436 
				false,  // 437 
				true,  // 438 
				true,  // 439 
				true,  // 440 
				true,  // 441 
				true,  // 442 
				true,  // 443 
				true,  // 444 
				true,  // 445 
				true,  // 446 
				true,  // 447 
				true,  // 448 
				true,  // 449 
				true,  // 450 
				true,  // 451 
				true,  // 452 
				true,  // 453 
				true,  // 454 
				true,  // 455 
				true,  // 456 
				true,  // 457 
				true,  // 458 
				true,  // 459 
				true,  // 460 
				true,  // 461 
				true,  // 462 
				true,  // 463 
				true,  // 464 
				true,  // 465 
				true,  // 466 
				true,  // 467 
				true,  // 468 
				true,  // 469 
				true,  // 470 
				true,  // 471 
				true,  // 472 
				true,  // 473 
				false,  // 474 
				false,  // 475 
				true,  // 476 
				true,  // 477 
				true,  // 478 
				true,  // 479 
				true,  // 480 
				true,  // 481 
				true,  // 482 
				true,  // 483 
				true,  // 484 
				true,  // 485 
				true,  // 486 
				true,  // 487 
				true,  // 488 
				true,  // 489 
				true,  // 490 
				true,  // 491 
				true,  // 492 
				true,  // 493 
				true,  // 494 
				true,  // 495 
				true,  // 496 
				true,  // 497 
				true,  // 498 
				true,  // 499 
				true,  // 500 
				true,  // 501 
				true,  // 502 
				true,  // 503 
				true,  // 504 
				true,  // 505 
				true,  // 506 
				true,  // 507 
				true,  // 508 
				true,  // 509 
				true,  // 510 
				true,  // 511 
				true,  // 512 
				true,  // 513 
				true,  // 514 
				true,  // 515 
				true,  // 516 
				true,  // 517 
				true,  // 518 
				true,  // 519 
				true,  // 520 
				true,  // 521 
				true,  // 522 
				true,  // 523 
				true,  // 524 
				true,  // 525 
				true,  // 526 
				true,  // 527 
				true,  // 528 
				true,  // 529 
				true,  // 530 
				true,  // 531 
				true,  // 532 
				true,  // 533 
				true,  // 534 
				true,  // 535 
				true,  // 536 
				true,  // 537 
				true,  // 538 
				true,  // 539 
				true,  // 540 
				true,  // 541 
				true,  // 542 
				true,  // 543 
				true,  // 544 
				false,  // 545 
				false,  // 546 
				false,  // 547 
				false,  // 548 
				false,  // 549 
				false,  // 550 
				false,  // 551 
				false,  // 552 
				false,  // 553 
				false,  // 554 
				false,  // 555 
				false,  // 556 
				true,  // 557 
				true,  // 558 
				true,  // 559 
				true,  // 560 
				true,  // 561 
				true,  // 562 
				true,  // 563 
				true,  // 564 
				true,  // 565 
				false,  // 566 
				false,  // 567 
				false,  // 568 
				true,  // 569 
				true,  // 570 
				true,  // 571 
				true,  // 572 
				true,  // 573 
				true,  // 574 
				true,  // 575 
				true,  // 576 
				true,  // 577 
				true,  // 578 
				true,  // 579 
				true,  // 580 
				true,  // 581 
				true,  // 582 
				true,  // 583 
				true,  // 584 
				true,  // 585 
				true,  // 586 
				true,  // 587 
				true,  // 588 
				true,  // 589 
				true,  // 590 
				true,  // 591 
				true,  // 592 
				true,  // 593 
				true,  // 594 
				true,  // 595 
				true,  // 596 
				true,  // 597 
				true,  // 598 
				true,  // 599 
				true,  // 600 
				true,  // 601 
				true,  // 602 
				true,  // 603 
				true,  // 604 
				true,  // 605 
				true,  // 606 
				true,  // 607 
				true,  // 608 
				true,  // 609 
				true,  // 610 
				true,  // 611 
				false,  // 612 
				false,  // 613 
				false,  // 614 
				false,  // 615 
				false,  // 616 
				false,  // 617 
				false,  // 618 
				false,  // 619 
				false,  // 620 
				false,  // 621 
				false,  // 622 
				false,  // 623 
				false,  // 624 
				false,  // 625 
				false,  // 626 
				false,  // 627 
				false,  // 628 
				false,  // 629 
				false,  // 630 
				false,  // 631 
				false,  // 632 
				false,  // 633 
				false };  // 634 
		
		Atom.resetNumberOfAtoms();
		AtomList newList = new AtomList();
		int templateCounter = 0;
		int queryCounter = 0;
		for (int resC=0 ; resC<templateSeq.length() ; resC++) {
			if (templateSeq.charAt(resC)!='-') {
				templateCounter++;
			}
			if (querySeq.charAt(resC)!='-') {
				queryCounter++;
			}
			if ((template.findAtomInList("CA", templateCounter)!=null) &
					(templateSeq.charAt(resC)!='-') &
					(querySeq.charAt(resC)!='-')) {
				Atom atom;
				// Doing N
				atom = template.findAtomInList("N", templateCounter);
				newList.add(new Atom(atom.x(), atom.y(), atom.z(),
						   "N", Residue.one2three(querySeq.charAt(resC)), queryCounter, -1));
				// Doing CA
				atom = template.findAtomInList("CA", templateCounter);
				newList.add(new Atom(atom.x(), atom.y(), atom.z(),
						   "CA", Residue.one2three(querySeq.charAt(resC)), queryCounter, -1));
				// Doing C
				atom = template.findAtomInList("C", templateCounter);
				newList.add(new Atom(atom.x(), atom.y(), atom.z(),
						   "C", Residue.one2three(querySeq.charAt(resC)), queryCounter, -1));
				// Doing O
				atom = template.findAtomInList("O", templateCounter);
				newList.add(new Atom(atom.x(), atom.y(), atom.z(),
						   "O", Residue.one2three(querySeq.charAt(resC)), queryCounter, -1));
			}
		}
		
		// Now putting a SCWRLED subunit
		AtomList newListSCWRLtmp = new AtomList("/Users/nirka/TRiC/New_Work_on_TRiC_14_8_2011/Gunnar/HM_Yeast_"+queryUnit+".scwrl.pdb");
		AtomList newListSCWRL = new AtomList();
		boolean[] resNumToTake = new boolean[residueInAllSeqs.length];
		for (int resC=0 ; resC<templateSeq.length() ; resC++) {
			resNumToTake[resC] = false;
		}
		queryCounter = 0;
		for (int resC=0 ; resC<templateSeq.length() ; resC++) {
			if (querySeq.charAt(resC)!='-') {
				queryCounter++;
				resNumToTake[queryCounter] = residueInAllSeqs[resC];
			}
		}
		for (int atomC=0 ; atomC<newListSCWRLtmp.size() ; atomC++) {
			if (resNumToTake[newListSCWRLtmp.atomAt(atomC).residueNumber()]) {
				newListSCWRL.add(newListSCWRLtmp.atomAt(atomC));
			}
		}
		AtomList equaQuery = filterDomainAccordingTo1Q3R(newListSCWRL, queryUnit , 'E');
		AtomList middleQuery = filterDomainAccordingTo1Q3R(newListSCWRL, queryUnit , 'M');
		AtomList apicQuery = filterDomainAccordingTo1Q3R(newListSCWRL, queryUnit , 'A');
		GDTcalculator.alignBySubset(newList, equaQuery, 0.75);
		GDTcalculator.alignBySubset(newList, middleQuery, 0.75);
		GDTcalculator.alignBySubset(newList, apicQuery, 0.75);
		return newListSCWRL;		
	}

	
	protected static AtomList turnWillisonTo1Q3R(AtomList oldList, String chainsTop, String genesTop, String chainsBot, String genesBot, boolean isE) {
		TricYeastAlignment alignments = new TricYeastAlignment();
		AtomList newList = new AtomList();
		// Doing top ring first
		for (int chain=0 ; chain<8 ; chain++) {
			AtomList chainAtoms = oldList.chainFilter(chainsTop.charAt(chain)+"");
			String alignment1Q3R = alignments.getAlignment("K");
			String alignmentChain = alignments.getAlignment(genesTop.charAt(chain)+"");
			int res1Q3RCounter=0;
			int resChainCounter=0;
			for (int res1Q3R=0; res1Q3R<alignment1Q3R.length() ; res1Q3R++) {
				if (alignment1Q3R.charAt(res1Q3R)!='-') {
					res1Q3RCounter++;
				}
				if (alignmentChain.charAt(res1Q3R)!='-') {
					resChainCounter++;
				}
				if ((alignment1Q3R.charAt(res1Q3R)!='-') & (alignmentChain.charAt(res1Q3R)!='-')) {
					Atom atom;
					if (isE) {
						atom = chainAtoms.findAtomInList("CA", resChainCounter+1000);
					}
					else {
						atom = chainAtoms.findAtomInList("CA", resChainCounter);
					}
					if (atom!=null) {
						Atom newAtom = new Atom(atom.x(), atom.y(), atom.z(),
								   "CA", Residue.one2three(alignment1Q3R.charAt(res1Q3R)), res1Q3RCounter, -1);
						newAtom.setChain(""+chainsTop.charAt(chain));
						newList.add(newAtom); 
					}
				}
			}					
		}
		// Doing bottom ring second
		for (int chain=0 ; chain<8 ; chain++) {
			AtomList chainAtoms = oldList.chainFilter(chainsBot.charAt(chain)+"");
			String alignment1Q3R = alignments.getAlignment("K");
			String alignmentChain = alignments.getAlignment(genesBot.charAt(chain)+"");
			int res1Q3RCounter=0;
			int resChainCounter=0;
			for (int res1Q3R=0; res1Q3R<alignment1Q3R.length() ; res1Q3R++) {
				if (alignment1Q3R.charAt(res1Q3R)!='-') {
					res1Q3RCounter++;
				}
				if (alignmentChain.charAt(res1Q3R)!='-') {
					resChainCounter++;
				}
				if ((alignment1Q3R.charAt(res1Q3R)!='-') & (alignmentChain.charAt(res1Q3R)!='-')) {
					Atom atom;
					if (isE) {
						atom = chainAtoms.findAtomInList("CA", resChainCounter+1000);
					}
					else {
						atom = chainAtoms.findAtomInList("CA", resChainCounter);
					}
					if (atom!=null) {
						Atom newAtom = new Atom(atom.x(), atom.y(), atom.z(),
								   "CA", Residue.one2three(alignment1Q3R.charAt(res1Q3R)), res1Q3RCounter, -1);
						newAtom.setChain(""+chainsBot.charAt(chain));
						newList.add(newAtom); 
					}
				}
			}					
		}
		return newList;
	}

	protected static AtomList filterDomainAccordingTo1Q3R(AtomList fullList, char subunit , char domain) {
		int[][] parsingQ3R; 
		TricYeastAlignment alignments = new TricYeastAlignment();
		switch (domain) {
    	case 'E': // Equatorial 
    		int[][] newParse = {{1,149} , {404,535}};
    		parsingQ3R = newParse;
    		break;
    	case 'M': // Middle 
    		int[][] newParse1 = {{150,217} , {369,403}};
    		parsingQ3R = newParse1;
    		break;
    	case 'A': // Apical 
    		int[][] newParse2 = {{218,368}};
    		parsingQ3R = newParse2;
    		break;
//    	case 'A': // Apical 
//    		int[][] newParse2 = {{218,247} , {278,368}};
//    		parsingQ3R = newParse2;
//    		break;
//    	case 'C': // Cap 
//    		int[][] newParse3 = {{248,277}};
//    		parsingQ3R = newParse3;
//    		break;
    	default:
    		throw new RuntimeException("Invalid domain letter {E,M,A}");
//    		throw new RuntimeException("Invalid domain letter {E,M,A,C}");
		}
		
		int[][] parsing = new int[parsingQ3R.length][2];
		for (int c=0 ; c<parsingQ3R.length ; c++) {
			parsing[c][0] = alignments.getNewResNum('K', parsingQ3R[c][0], subunit);
			parsing[c][1] = alignments.getNewResNum('K', parsingQ3R[c][1], subunit);
		}		
		return filterDomain(fullList, parsing);
	}
	
	
	protected static AtomList filterDomain(AtomList fullList, int[][] parsing) {
		AtomList newList = new AtomList();
		for (int c=0 ; c<fullList.size() ; c++) {
			int resNum = fullList.atomAt(c).residueNumber();
			for (int segID=0 ; segID<parsing.length ; segID++ ) {
				if ((resNum>=parsing[segID][0]) & (resNum<=parsing[segID][1])) {
					newList.add(fullList.atomAt(c));
				}
			}
		}
		return newList;
	}
	
	
	protected static void init(String[] args) {
		int zvl = ALA; // force the reading of "meshi.parameters.Residues"
		zvl = ACA;// force the reading of "meshi.parameters.AtomTypes"
		initRandom(333);
	}	

}
