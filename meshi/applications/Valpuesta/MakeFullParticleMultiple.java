package meshi.applications.Valpuesta;

import java.io.BufferedWriter;
import java.io.FileWriter;

import meshi.applications.TriC.TricAlignment;
import meshi.applications.TriC.TricYeastAlignment;
import meshi.applications.prediction.GDTcalculator;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Residue;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.MeshiProgram;
import meshi.util.file.File2StringArray;

public class MakeFullParticleMultiple extends MeshiProgram implements Residues,AtomTypes {
	
	public static void main(String[] args) throws Exception {
		init(args);
		String[] allModels = File2StringArray.f2a("Top_100_for_Gunnar_models.txt");
		for (int c=0; c<allModels.length ; c++) {
			makeSingle(allModels[c].substring(0, 8), allModels[c].substring(8, 16), c+1, allModels[c].substring(17, 25));
		}
	}
		
	public static void makeSingle(String topTrue, String botTrue, int modelNumber, String Rstring) throws Exception {
		// Defining the arrangement 
//		String topTrue   = "ZQHEBDAG"; // OMS
//		String topTrue   = "ZEAHDQGB"; // Willison
//		String topTrue   = "AZBGQDEH"; // Yao
		String chainsTop = "ABCDEFGH";
		
//		String botTrue = "ZQHEBDAG"; // OMS
//		String botTrue = "ZEAHDQGB"; // Willison
//		String botTrue = "AZBGQDEH"; // Yao
		String chainsBot = "PONMLKJI";  // P chain is right under A chain

		AtomList willison1Q3R = new AtomList("2XSM_as_1Q3R.pdb");
		Atom.resetNumberOfAtoms();
		TricAlignment alignments = new TricAlignment();

		
		for (int delta1=0 ; delta1<1 ; delta1++) { // Just 1 structure
			String outputFile;
			if (modelNumber>99) {
				outputFile = "Top_Models_for_Gunnar/model_"+modelNumber+"_"+topTrue+"_"+botTrue+"_"+Rstring+".pdb";
			}
			else if (modelNumber>9) {
				outputFile = "Top_Models_for_Gunnar/model_0"+modelNumber+"_"+topTrue+"_"+botTrue+"_"+Rstring+".pdb";
			}
			else {
				outputFile = "Top_Models_for_Gunnar/model_00"+modelNumber+"_"+topTrue+"_"+botTrue+"_"+Rstring+".pdb";
			}
			BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile));
			bw.write("CRYST1  272.700  313.500  158.300  90.00  90.00  90.00 P 21 21 2     8          \n");                         
			modelNumber++;
			// Doing top ring
			for (int pos1=0 ; pos1<8 ; pos1++) { // Position on Willison's PDB
				char chainID = chainsTop.charAt(pos1);
				char toPut = topTrue.charAt((pos1+8-delta1) % 8);
				AtomList template = willison1Q3R.chainFilter(chainID+"");
				Atom.resetNumberOfAtoms();
				System.out.println("\n\n" + toPut + " in chain " + chainID + "  ::: " + template.size());
				AtomList finalUnit = HMunit(alignments.getAlignment(toPut+""), toPut,alignments.getAlignment("K"), 'K',template);
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
				AtomList template = willison1Q3R.chainFilter(chainID+"");
				Atom.resetNumberOfAtoms();
				System.out.println("\n\n" + toPut + " in chain " + chainID + "  ::: " + template.size());
				AtomList finalUnit = HMunit(alignments.getAlignment(toPut+""), toPut,alignments.getAlignment("K"), 'K',template);
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

	
	protected static AtomList HMunit(String querySeq,char queryUnit , String templateSeq, char templateUnit, AtomList template) {
		boolean[] residueInAllSeqs = {false,  // -   
				false,  // - 1 
				false,  // - 2 
				false,  // - 3 
				false,  // - 4 
				false,  // - 5 
				false,  // - 6 
				false,  // - 7 
				false,  // - 8 
				false,  // - 9 
				false,  // - 10 
				false,  // - 11 
				false,  // - 12 
				false,  // - 13 
				false,  // - 14 
				false,  // - 15 
				false,  // - 16 
				false,  // - 17 
				false,  // - 18 
				false,  // - 19 
				false,  // - 20 
				false,  // - 21 
				false,  // M 22 
				false,  // A 23 
				false,  // Q 24 
				false,  // L 25 
				false,  // S 26 
				false,  // G 27 
				false,  // Q 28 
				true,  // P 29 
				true,  // V 30 
				true,  // V 31 
				true,  // I 32 
				true,  // L 33 
				true,  // P 34 
				true,  // E 35 
				true,  // G 36 
				true,  // T 37 
				true,  // Q 38 
				true,  // R 39 
				true,  // Y 40 
				true,  // V 41 
				true,  // G 42 
				true,  // R 43 
				true,  // D 44 
				true,  // A 45 
				true,  // Q 46 
				true,  // R 47 
				true,  // L 48 
				true,  // N 49 
				true,  // I 50 
				true,  // L 51 
				true,  // A 52 
				true,  // A 53 
				true,  // R 54 
				true,  // I 55 
				true,  // I 56 
				true,  // A 57 
				true,  // E 58 
				true,  // T 59 
				true,  // V 60 
				true,  // R 61 
				true,  // T 62 
				true,  // T 63 
				true,  // L 64 
				true,  // G 65 
				true,  // P 66 
				true,  // K 67 
				true,  // G 68 
				true,  // M 69 
				true,  // D 70 
				true,  // K 71 
				true,  // M 72 
				true,  // L 73 
				true,  // V 74 
				true,  // D 75 
				true,  // S 76 
				false,  // - 77 
				false,  // - 78 
				true,  // L 79 
				true,  // G 80 
				true,  // D 81 
				true,  // I 82 
				true,  // V 83 
				true,  // V 84 
				true,  // T 85 
				true,  // N 86 
				true,  // D 87 
				true,  // G 88 
				true,  // A 89 
				true,  // T 90 
				true,  // I 91 
				true,  // L 92 
				true,  // D 93 
				true,  // K 94 
				true,  // I 95 
				true,  // D 96 
				true,  // L 97 
				true,  // Q 98 
				true,  // H 99 
				true,  // P 100 
				true,  // A 101 
				true,  // A 102 
				true,  // K 103 
				true,  // M 104 
				true,  // M 105 
				true,  // V 106 
				true,  // E 107 
				true,  // V 108 
				true,  // A 109 
				true,  // K 110 
				true,  // T 111 
				true,  // Q 112 
				true,  // D 113 
				true,  // K 114 
				true,  // E 115 
				true,  // A 116 
				true,  // G 117 
				true,  // D 118 
				true,  // G 119 
				true,  // T 120 
				true,  // T 121 
				true,  // T 122 
				true,  // A 123 
				true,  // V 124 
				true,  // V 125 
				true,  // I 126 
				true,  // A 127 
				true,  // G 128 
				true,  // E 129 
				true,  // L 130 
				true,  // L 131 
				true,  // R 132 
				true,  // K 133 
				true,  // A 134 
				true,  // E 135 
				true,  // E 136 
				true,  // L 137 
				true,  // L 138 
				false,  // - 139 
				true,  // D 140 
				true,  // Q 141 
				true,  // N 142 
				true,  // I 143 
				true,  // H 144 
				true,  // P 145 
				true,  // S 146 
				true,  // I 147 
				true,  // I 148 
				true,  // T 149 
				true,  // K 150 
				true,  // G 151 
				true,  // Y 152 
				true,  // A 153 
				true,  // L 154 
				true,  // A 155 
				true,  // A 156 
				true,  // E 157 
				true,  // K 158 
				true,  // A 159 
				true,  // Q 160 
				true,  // E 161 
				true,  // I 162 
				true,  // L 163 
				true,  // D 164 
				true,  // E 165 
				true,  // I 166 
				true,  // A 167 
				true,  // I 168 
				true,  // R 169 
				true,  // V 170 
				true,  // D 171 
				false,  // - 172 
				false,  // - 173 
				false,  // - 174 
				false,  // - 175 
				false,  // P 176 
				true,  // D 177 
				true,  // D 178 
				true,  // E 179 
				true,  // E 180 
				true,  // T 181 
				true,  // L 182 
				true,  // L 183 
				true,  // K 184 
				true,  // I 185 
				true,  // A 186 
				true,  // A 187 
				true,  // T 188 
				true,  // S 189 
				true,  // I 190 
				true,  // T 191 
				true,  // G 192 
				true,  // K 193 
				true,  // N 194 
				true,  // A 195 
				true,  // E 196 
				true,  // S 197 
				true,  // H 198 
				true,  // K 199 
				true,  // E 200 
				true,  // L 201 
				true,  // L 202 
				true,  // A 203 
				true,  // K 204 
				true,  // L 205 
				true,  // A 206 
				true,  // V 207 
				true,  // E 208 
				true,  // A 209 
				true,  // V 210 
				true,  // K 211 
				true,  // Q 212 
				true,  // V 213 
				true,  // A 214 
				true,  // E 215 
				false,  // K 216 
				false,  // K 217 
				false,  // D 218 
				false,  // - 219 
				false,  // - 220 
				false,  // - 221 
				false,  // - 222 
				false,  // - 223 
				false,  // - 224 
				false,  // - 225 
				false,  // G 226 
				false,  // K 227 
				false,  // Y 228 
				true,  // V 229 
				true,  // V 230 
				true,  // D 231 
				true,  // L 232 
				true,  // D 233 
				true,  // N 234 
				true,  // I 235 
				true,  // K 236 
				true,  // F 237 
				true,  // E 238 
				true,  // K 239 
				true,  // K 240 
				true,  // A 241 
				true,  // G 242 
				true,  // E 243 
				true,  // G 244 
				true,  // V 245 
				true,  // E 246 
				true,  // E 247 
				true,  // S 248 
				true,  // E 249 
				true,  // L 250 
				true,  // V 251 
				true,  // R 252 
				true,  // G 253 
				true,  // V 254 
				true,  // V 255 
				true,  // I 256 
				true,  // D 257 
				false,  // - 258 
				false,  // - 259 
				false,  // - 260 
				false,  // K 261 
				false,  // E 262 
				true,  // V 263 
				true,  // V 264 
				true,  // H 265 
				true,  // P 266 
				true,  // R 267 
				false,  // - 268 
				false,  // - 269 
				true,  // M 270 
				true,  // P 271 
				true,  // K 272 
				true,  // R 273 
				true,  // V 274 
				true,  // E 275 
				false,  // - 276 
				false,  // - 277 
				false,  // - 278 
				false,  // - 279 
				false,  // - 280 
				false,  // - 281 
				true,  // N 282 
				true,  // A 283 
				true,  // K 284 
				true,  // I 285 
				true,  // A 286 
				true,  // L 287 
				true,  // I 288 
				true,  // N 289 
				true,  // E 290 
				true,  // A 291 
				true,  // L 292 
				true,  // E 293 
				true,  // V 294 
				true,  // K 295 
				true,  // K 296 
				true,  // T 297 
				true,  // E 298 
				false,  // - 299 
				true,  // T 300 
				true,  // D 301 
				true,  // A 302 
				true,  // K 303 
				true,  // I 304 
				true,  // N 305 
				true,  // I 306 
				true,  // T 307 
				true,  // S 308 
				true,  // P 309 
				true,  // D 310 
				true,  // Q 311 
				true,  // L 312 
				true,  // M 313 
				true,  // S 314 
				true,  // F 315 
				true,  // L 316 
				true,  // E 317 
				true,  // Q 318 
				true,  // E 319 
				true,  // E 320 
				true,  // K 321 
				true,  // M 322 
				true,  // L 323 
				true,  // K 324 
				true,  // D 325 
				true,  // M 326 
				true,  // V 327 
				true,  // D 328 
				true,  // H 329 
				true,  // I 330 
				true,  // A 331 
				true,  // Q 332 
				true,  // T 333 
				true,  // G 334 
				false,  // - 335 
				false,  // - 336 
				false,  // - 337 
				false,  // - 338 
				false,  // - 339 
				false,  // - 340 
				false,  // - 341 
				false,  // - 342 
				false,  // - 343 
				false,  // - 344 
				true,  // A 345 
				true,  // N 346 
				true,  // V 347 
				true,  // V 348 
				true,  // F 349 
				true,  // V 350 
				true,  // Q 351 
				false,  // - 352 
				false,  // - 353 
				false,  // - 354 
				false,  // - 355 
				false,  // - 356 
				true,  // K 357 
				true,  // G 358 
				true,  // I 359 
				true,  // D 360 
				true,  // D 361 
				true,  // L 362 
				true,  // A 363 
				true,  // Q 364 
				true,  // H 365 
				true,  // Y 366 
				true,  // L 367 
				true,  // A 368 
				true,  // K 369 
				true,  // Y 370 
				true,  // G 371 
				true,  // I 372 
				true,  // M 373 
				true,  // A 374 
				true,  // V 375 
				true,  // R 376 
				true,  // R 377 
				true,  // V 378 
				true,  // K 379 
				true,  // K 380 
				true,  // S 381 
				true,  // D 382 
				true,  // M 383 
				true,  // E 384 
				true,  // K 385 
				true,  // L 386 
				true,  // A 387 
				true,  // K 388 
				true,  // A 389 
				true,  // T 390 
				true,  // G 391 
				true,  // A 392 
				true,  // K 393 
				true,  // I 394 
				true,  // V 395 
				true,  // T 396 
				true,  // N 397 
				true,  // V 398 
				true,  // K 399 
				true,  // D 400 
				true,  // L 401 
				true,  // T 402 
				false,  // - 403 
				false,  // - 404 
				false,  // - 405 
				false,  // - 406 
				false,  // - 407 
				false,  // - 408 
				true,  // P 409 
				true,  // E 410 
				true,  // D 411 
				true,  // L 412 
				true,  // G 413 
				true,  // Y 414 
				true,  // A 415 
				true,  // E 416 
				false,  // - 417 
				true,  // V 418 
				true,  // V 419 
				true,  // E 420 
				true,  // E 421 
				true,  // R 422 
				true,  // K 423 
				true,  // L 424 
				true,  // A 425 
				false,  // - 426 
				false,  // - 427 
				true,  // G 428 
				true,  // E 429 
				true,  // N 430 
				true,  // M 431 
				true,  // I 432 
				true,  // F 433 
				true,  // V 434 
				true,  // E 435 
				true,  // G 436 
				true,  // C 437 
				true,  // K 438 
				true,  // N 439 
				false,  // - 440 
				false,  // - 441 
				false,  // - 442 
				true,  // P 443 
				true,  // K 444 
				true,  // A 445 
				true,  // V 446 
				true,  // T 447 
				true,  // I 448 
				true,  // L 449 
				true,  // I 450 
				true,  // R 451 
				true,  // G 452 
				true,  // G 453 
				true,  // T 454 
				true,  // E 455 
				true,  // H 456 
				true,  // V 457 
				true,  // I 458 
				true,  // D 459 
				true,  // E 460 
				true,  // V 461 
				true,  // E 462 
				true,  // R 463 
				true,  // A 464 
				true,  // L 465 
				true,  // E 466 
				true,  // D 467 
				true,  // A 468 
				true,  // V 469 
				true,  // K 470 
				true,  // V 471 
				true,  // V 472 
				true,  // K 473 
				true,  // D 474 
				true,  // V 475 
				true,  // M 476 
				true,  // E 477 
				true,  // D 478 
				false,  // - 479 
				false,  // - 480 
				true,  // G 481 
				true,  // A 482 
				true,  // V 483 
				true,  // L 484 
				true,  // P 485 
				true,  // A 486 
				true,  // G 487 
				true,  // G 488 
				true,  // A 489 
				true,  // P 490 
				true,  // E 491 
				true,  // I 492 
				true,  // E 493 
				true,  // L 494 
				true,  // A 495 
				true,  // I 496 
				true,  // R 497 
				true,  // L 498 
				true,  // D 499 
				true,  // E 500 
				true,  // Y 501 
				true,  // A 502 
				true,  // K 503 
				true,  // Q 504 
				true,  // V 505 
				true,  // G 506 
				false,  // - 507 
				false,  // - 508 
				false,  // - 509 
				true,  // G 510 
				true,  // K 511 
				true,  // E 512 
				true,  // A 513 
				true,  // L 514 
				true,  // A 515 
				true,  // I 516 
				true,  // E 517 
				true,  // N 518 
				true,  // F 519 
				true,  // A 520 
				true,  // D 521 
				true,  // A 522 
				true,  // L 523 
				true,  // K 524 
				true,  // I 525 
				true,  // I 526 
				true,  // P 527 
				true,  // K 528 
				true,  // T 529 
				true,  // L 530 
				true,  // A 531 
				true,  // E 532 
				true,  // N 533 
				true,  // A 534 
				true,  // G 535 
				true,  // L 536 
				true,  // D 537 
				true,  // T 538 
				true,  // V 539 
				true,  // E 540 
				true,  // M 541 
				true,  // L 542 
				true,  // V 543 
				true,  // K 544 
				true,  // V 545 
				true,  // I 546 
				true,  // S 547 
				true,  // E 548 
				true,  // H 549 
				true,  // K 550 
				true,  // N 551 
				true,  // R 552 
				false,  // - 553 
				false,  // - 554 
				false,  // - 555 
				false,  // - 556 
				false,  // - 557 
				false,  // - 558 
				false,  // - 559 
				false,  // - 560 
				false,  // - 561 
				false,  // - 562 
				false,  // - 563 
				false,  // - 564 
				true,  // G 565 
				true,  // L 566 
				true,  // G 567 
				true,  // I 568 
				true,  // G 569 
				true,  // I 570 
				true,  // D 571 
				true,  // V 572 
				true,  // F 573 
				false,  // - 574 
				false,  // - 575 
				false,  // - 576 
				true,  // E 577 
				true,  // G 578 
				true,  // K 579 
				true,  // P 580 
				true,  // A 581 
				true,  // D 582 
				true,  // M 583 
				true,  // L 584 
				true,  // E 585 
				true,  // K 586 
				true,  // G 587 
				true,  // I 588 
				true,  // I 589 
				true,  // E 590 
				true,  // P 591 
				true,  // L 592 
				true,  // R 593 
				true,  // V 594 
				true,  // K 595 
				true,  // K 596 
				true,  // Q 597 
				true,  // A 598 
				true,  // I 599 
				true,  // K 600 
				true,  // S 601 
				true,  // A 602 
				true,  // S 603 
				true,  // E 604 
				true,  // A 605 
				true,  // A 606 
				true,  // I 607 
				true,  // M 608 
				true,  // I 609 
				true,  // L 610 
				true,  // R 611 
				true,  // I 612 
				true,  // D 613 
				true,  // D 614 
				true,  // V 615 
				true,  // I 616 
				true,  // A 617 
				true,  // A 618 
				true,  // K 619 
				false,  // A 620 
				false,  // T 621 
				false,  // K 622 
				false,  // P 623 
				false,  // E 624 
				false,  // G 625 
				false,  // G 626 
				false,  // Q 627 
				false,  // G 628 
				false,  // G 629 
				false,  // G 630 
				false,  // M 631 
				false,  // P 632 
				false,  // G 633 
				false,  // G 634 
				false,  // M 635 
				false,  // G 636 
				false,  // G 637 
				false,  // M 638 
				false,  // G 639 
				false,  // M 640 
				false,  // G 641 
				false};  // M 642 
		
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
				// Doing CA
				atom = template.findAtomInList("CA", templateCounter);
				newList.add(new Atom(atom.x(), atom.y(), atom.z(),
						   "CA", Residue.one2three(querySeq.charAt(resC)), queryCounter, -1));
			}
			else if (template.findAtomInList("CA", templateCounter)==null) {
				residueInAllSeqs[resC] = false;
			}
		}
		
		// Now putting a SCWRLED subunit
		AtomList newListSCWRLtmp = new AtomList("C:\\Users\\Nir\\TRiC\\Valpu\\CCT"+queryUnit+"_BOV_scwrl.pdb");
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
//		System.out.println("Sizes: " + equaQuery.CAFilter().size() + " " + middleQuery.CAFilter().size() + " " + apicQuery.CAFilter().size() + "\n"); 
//		if ((middleQuery.CAFilter().size()<11)&(middleQuery.CAFilter().size()>0)) {
//			middleQuery.print();
//			System.out.println("\n\n" + queryUnit);
//			System.exit(0);
//		}
//		if ((apicQuery.CAFilter().size()<11)&(apicQuery.CAFilter().size()>0)) {
//			apicQuery.print();
//			System.out.println("\n\n" + queryUnit);
//			System.exit(0);
//		}
		if (equaQuery.size()>0 )
			GDTcalculator.alignBySubset(newList, equaQuery, 0.75);
		if (middleQuery.size()>0 )
			GDTcalculator.alignBySubset(newList, middleQuery, 0.75);
		if (apicQuery.size()>0 )
			GDTcalculator.alignBySubset(newList, apicQuery, 0.75);
		return newListSCWRL;		
	}

	

	protected static AtomList filterDomainAccordingTo1Q3R(AtomList fullList, char subunit , char domain) {
		int[][] parsingQ3R; 
		TricAlignment alignments = new TricAlignment();
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


 