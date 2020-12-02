package uniandes.algorithms.tr;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import ngsep.sequences.SimpleEditDistanceMeasure;
import ngsep.sequences.PairwiseAlignmentAffineGap;

public class TRFCandidateSelector {

	private String sequence;

	// probe p, list of occurrences
	private HashMap<String, ArrayList<Integer>> historyList;

	// distance d, list of matches (last index)
	private HashMap<Integer, ArrayList<Integer>> distanceList;

	// probability of a match, can be 0.8 or 0.75
	private double matchProbability;

	// probability of indel, typically 0.1
	private double indelProbability;

	// list of candidates in case multiple tuple sizes are used
	private ArrayList<TandemRepeat> allCandidates;

	// max pattern size to look for, the bigger it is, the more time the
	// algorithm takes
	private int maxPatternSize;

	// minimum alignment score for a tandem repeat to be reported
	private int minimumAlignmentScore;

	// sequence name for this instance
	private String sequenceName;

	private int matchScore;

	private int missScore;

	// sum of data criteria for pm = 0.8, retrieved from
	// https://github.com/Benson-Genomics-Lab/TRF/blob/master/src/tr30dat.c
	int sumdata80[] = { 0, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 7, 8, 8, 9, 9, 6, 6,
			7, 7, 7, 8, 8, 9, 9, 10, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 16, 17, 17, 18, 18, 19, 19,
			20, 20, 21, 21, 22, 22, 23, 23, 24, 24, 24, 25, 25, 26, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 32, 32,
			33, 33, 34, 34, 35, 35, 36, 36, 37, 37, 38, 38, 39, 39, 40, 40, 41, 41, 42, 42, 43, 43, 44, 44, 45, 45, 46,
			46, 47, 47, 48, 48, 49, 49, 50, 50, 51, 51, 52, 52, 53, 53, 54, 54, 55, 55, 56, 56, 57, 57, 58, 58, 59, 59,
			60, 60, 61, 61, 62, 63, 63, 64, 64, 65, 65, 66, 66, 67, 67, 68, 68, 69, 43, 43, 43, 44, 44, 44, 45, 45, 46,
			46, 46, 47, 47, 47, 48, 48, 49, 49, 49, 50, 50, 50, 51, 51, 52, 52, 52, 53, 53, 53, 54, 54, 55, 55, 55, 56,
			56, 56, 57, 57, 58, 58, 58, 59, 59, 59, 60, 60, 61, 61, 61, 62, 62, 63, 63, 63, 64, 64, 64, 65, 65, 66, 66,
			66, 67, 67, 68, 68, 68, 69, 69, 69, 70, 70, 71, 71, 71, 72, 72, 73, 73, 73, 74, 74, 74, 75, 75, 76, 76, 76,
			77, 77, 78, 78, 78, 79, 79, 80, 80, 80, 81, 81, 81, 82, 82, 83, 83, 83, 84, 84, 85, 85, 85, 86, 86, 87, 87,
			87, 88, 88, 89, 89, 89, 90, 90, 90, 91, 91, 92, 92, 92, 93, 93, 94, 94, 94, 95, 95, 96, 96, 96, 97, 97, 98,
			98, 98, 99, 99, 100, 100, 100, 101, 101, 102, 102, 102, 103, 103, 103, 104, 104, 105, 105, 105, 106, 106,
			107, 107, 107, 108, 108, 109, 109, 109, 110, 110, 111, 111, 111, 112, 112, 113, 113, 113, 114, 114, 115,
			115, 115, 116, 116, 117, 117, 117, 118, 118, 119, 119, 119, 120, 120, 121, 121, 121, 122, 122, 123, 123,
			123, 124, 124, 125, 125, 125, 126, 126, 127, 127, 127, 128, 128, 129, 129, 129, 130, 130, 131, 131, 131,
			132, 132, 133, 133, 134, 134, 134, 135, 135, 136, 136, 136, 137, 137, 138, 138, 138, 139, 139, 140, 140,
			140, 141, 141, 142, 142, 142, 143, 143, 144, 144, 144, 145, 145, 146, 146, 146, 147, 147, 148, 148, 148,
			149, 149, 150, 150, 151, 151, 151, 152, 152, 153, 153, 153, 154, 154, 155, 155, 155, 156, 156, 157, 157,
			157, 158, 158, 159, 159, 159, 160, 160, 161, 161, 161, 162, 162, 163, 163, 164, 164, 164, 165, 165, 166,
			166, 166, 167, 167, 168, 168, 168, 169, 169, 170, 170, 170, 171, 171, 172, 172, 173, 173, 173, 174, 174,
			175, 175, 175, 176, 176, 177, 177, 177, 178, 178, 179, 179, 179, 180, 180, 181, 181, 182, 182, 182, 183,
			183, 184, 184, 184, 185, 185, 186, 186, 186, 187, 187, 188, 188, 189, 189, 189, 190, 190, 191, 191, 191,
			192, 192, 193, 193, 193, 194, 194, 195, 195, 196, 196, 196, 197, 197, 198, 198, 198, 199, 199, 200, 200,
			200, 201, 201, 202, 202, 203, 203, 203, 204, 204, 205, 205, 205, 206, 206, 207, 207, 207, 208, 208, 209,
			209, 210, 210, 210, 211, 211, 212, 212, 212, 213, 213, 214, 214, 215, 215, 215, 216, 216, 217, 217, 217,
			218, 218, 219, 219, 219, 220, 220, 221, 221, 222, 222, 222, 223, 223, 224, 224, 224, 225, 225, 226, 226,
			227, 227, 227, 228, 228, 229, 229, 229, 230, 230, 231, 231, 232, 232, 232, 233, 233, 234, 234, 234, 235,
			235, 236, 236, 237, 237, 237, 238, 238, 239, 239, 239, 240, 240, 241, 241, 242, 242, 242, 243, 243, 244,
			244, 244, 245, 245, 246, 246, 246, 247, 247, 248, 248, 249, 249, 249, 250, 250, 251, 251, 252, 252, 252,
			253, 253, 254, 254, 254, 255, 255, 256, 256, 257, 257, 257, 258, 258, 259, 259, 259, 260, 260, 261, 261,
			262, 262, 262, 263, 263, 264, 264, 264, 265, 265, 266, 266, 267, 267, 267, 268, 268, 269, 269, 269, 270,
			270, 271, 271, 272, 272, 272, 273, 273, 274, 274, 274, 275, 275, 276, 276, 277, 277, 277, 278, 278, 279,
			279, 280, 280, 280, 281, 281, 282, 282, 282, 283, 283, 284, 284, 285, 285, 285, 286, 286, 287, 287, 287,
			288, 288, 289, 289, 290, 290, 290, 291, 291, 292, 292, 293, 293, 293, 294, 294, 295, 295, 295, 296, 296,
			297, 297, 298, 298, 298, 299, 299, 300, 300, 301, 301, 301, 302, 302, 303, 303, 303, 304, 304, 305, 305,
			306, 306, 306, 307, 307, 308, 308, 309, 309, 309, 310, 310, 311, 311, 311, 312, 312, 313, 313, 314, 314,
			314, 315, 315, 316, 316, 317, 317, 317, 318, 318, 319, 319, 319, 320, 320, 321, 321, 322, 322, 322, 323,
			323, 324, 324, 325, 325, 325, 326, 326, 327, 327, 327, 328, 328, 329, 329, 330, 330, 330, 331, 331, 332,
			332, 333, 333, 333, 334, 334, 335, 335, 336, 336, 336, 337, 337, 338, 338, 338, 339, 339, 340, 340, 341,
			341, 341, 342, 342, 343, 343, 344, 344, 344, 345, 345, 346, 346, 347, 347, 347, 348, 348, 349, 349, 349,
			350, 350, 351, 351, 352, 352, 352, 353, 353, 354, 354, 355, 355, 355, 356, 356, 357, 357, 358, 358, 358,
			359, 359, 360, 360, 360, 361, 361, 362, 362, 363, 363, 363, 364, 364, 365, 365, 366, 366, 366, 367, 367,
			368, 368, 369, 369, 369, 370, 370, 371, 371, 372, 372, 372, 373, 373, 374, 374, 374, 375, 375, 376, 376,
			377, 377, 377, 378, 378, 379, 379, 380, 380, 380, 381, 381, 382, 382, 383, 383, 383, 384, 384, 385, 385,
			386, 386, 386, 387, 387, 388, 388, 388, 389, 389, 390, 390, 391, 391, 391, 392, 392, 393, 393, 394, 394,
			394, 395, 395, 396, 396, 397, 397, 397, 398, 398, 399, 399, 400, 400, 400, 401, 401, 402, 402, 403, 403,
			403, 404, 404, 405, 405, 406, 406, 406, 407, 407, 408, 408, 408, 409, 409, 410, 410, 411, 411, 411, 412,
			412, 413, 413, 414, 414, 414, 415, 415, 416, 416, 417, 417, 417, 418, 418, 419, 419, 420, 420, 420, 421,
			421, 422, 422, 423, 423, 423, 424, 424, 425, 425, 426, 426, 426, 427, 427, 428, 428, 429, 429, 429, 430,
			430, 431, 431, 432, 432, 432, 433, 433, 434, 434, 434, 435, 435, 436, 436, 437, 437, 437, 438, 438, 439,
			439, 440, 440, 440, 441, 441, 442, 442, 443, 443, 443, 444, 444, 445, 445, 446, 446, 446, 447, 447, 448,
			448, 449, 449, 449, 450, 450, 451, 451, 452, 452, 452, 453, 453, 454, 454, 455, 455, 455, 456, 456, 457,
			457, 458, 458, 458, 459, 459, 460, 460, 461, 461, 461, 462, 462, 463, 463, 464, 464, 464, 465, 465, 466,
			466, 467, 467, 467, 468, 468, 469, 469, 470, 470, 470, 471, 471, 472, 472, 473, 473, 473, 474, 474, 475,
			475, 476, 476, 476, 477, 477, 478, 478, 479, 479, 479, 480, 480, 481, 481, 482, 482, 482, 483, 483, 484,
			484, 485, 485, 485, 486, 486, 487, 487, 488, 488, 488, 489, 489, 490, 490, 491, 491, 491, 492, 492, 493,
			493, 494, 494, 494, 495, 495, 496, 496, 497, 497, 497, 498, 498, 499, 499, 500, 500, 500, 501, 501, 502,
			502, 503, 503, 503, 504, 504, 505, 505, 506, 506, 506, 507, 507, 508, 508, 509, 509, 509, 510, 510, 511,
			511, 512, 512, 512, 513, 513, 514, 514, 515, 515, 515, 516, 516, 517, 517, 518, 518, 518, 519, 519, 520,
			520, 521, 521, 521, 522, 522, 523, 523, 524, 524, 524, 525, 525, 526, 526, 527, 527, 527, 528, 528, 529,
			529, 530, 530, 530, 531, 531, 532, 532, 533, 533, 533, 534, 534, 535, 535, 536, 536, 536, 537, 537, 538,
			538, 539, 539, 539, 540, 540, 541, 541, 542, 542, 543, 543, 543, 544, 544, 545, 545, 546, 546, 546, 547,
			547, 548, 548, 549, 549, 549, 550, 550, 551, 551, 552, 552, 552, 553, 553, 554, 554, 555, 555, 555, 556,
			556, 557, 557, 558, 558, 558, 559, 559, 560, 560, 561, 561, 561, 562, 562, 563, 563, 564, 564, 564, 565,
			565, 566, 566, 567, 567, 567, 568, 568, 569, 569, 570, 570, 570, 571, 571, 572, 572, 573, 573, 574, 574,
			574, 575, 575, 576, 576, 577, 577, 577, 578, 578, 579, 579, 580, 580, 580, 581, 581, 582, 582, 583, 583,
			583, 584, 584, 585, 585, 586, 586, 586, 587, 587, 588, 588, 589, 589, 589, 590, 590, 591, 591, 592, 592,
			592, 593, 593, 594, 594, 595, 595, 595, 596, 596, 597, 597, 598, 598, 599, 599, 599, 600, 600, 601, 601,
			602, 602, 602, 603, 603, 604, 604, 605, 605, 605, 606, 606, 607, 607, 608, 608, 608, 609, 609, 610, 610,
			611, 611, 611, 612, 612, 613, 613, 614, 614, 614, 615, 615, 616, 616, 617, 617, 618, 618, 618, 619, 619,
			620, 620, 621, 621, 621, 622, 622, 623, 623, 624, 624, 624, 625, 625, 626, 626, 627, 627, 627, 628, 628,
			629, 629, 630, 630, 630, 631, 631, 632, 632, 633, 633, 634, 634, 634, 635, 635, 636, 636, 637, 637, 637,
			638, 638, 639, 639, 640, 640, 640, 641, 641, 642, 642, 643, 643, 643, 644, 644, 645, 645, 646, 646, 646,
			647, 647, 648, 648, 649, 649, 650, 650, 650, 651, 651, 652, 652, 653, 653, 653, 654, 654, 655, 655, 656,
			656, 656, 657, 657, 658, 658, 659, 659, 659, 660, 660, 661, 661, 662, 662, 662, 663, 663, 664, 664, 665,
			665, 666, 666, 666, 667, 667, 668, 668, 669, 669, 669, 670, 670, 671, 671, 672, 672, 672, 673, 673, 674,
			674, 675, 675, 675, 676, 676, 677, 677, 678, 678, 679, 679, 679, 680, 680, 681, 681, 682, 682, 682, 683,
			683, 684, 684, 685, 685, 685, 686, 686, 687, 687, 688, 688, 688, 689, 689, 690, 690, 691, 691, 692, 692,
			692, 693, 693, 694, 694, 695, 695, 695, 696, 696, 697, 697, 698, 698, 698, 699, 699, 700, 700, 701, 701,
			701, 702, 702, 703, 703, 704, 704, 705, 705, 705, 706, 706, 707, 707, 708, 708, 708, 709, 709, 710, 710,
			711, 711, 711, 712, 712, 713, 713, 714, 714, 715, 715, 715, 716, 716, 717, 717, 718, 718, 718, 719, 719,
			720, 720, 721, 721, 721, 722, 722, 723, 723, 724, 724, 724, 725, 725, 726, 726, 727, 727, 728, 728, 728,
			729, 729, 730, 730, 731, 731, 731, 732, 732, 733, 733, 734, 734, 734, 735, 735, 736, 736, 737, 737, 738,
			738, 738, 739, 739, 740, 740, 741, 741, 741, 742, 742, 743, 743, 744, 744, 744, 745, 745, 746, 746, 747,
			747, 748, 748, 748, 749, 749, 750, 750, 751, 751, 751, 752, 752, 753, 753, 754, 754, 754, 755, 755, 756,
			756, 757, 757, 758, 758, 758, 759, 759, 760, 760, 761, 761, 761, 762, 762, 763, 763, 764, 764, 764, 765,
			765, 766, 766, 767, 767, 768, 768, 768, 769, 769, 770, 770, 771, 771, 771, 772, 772, 773, 773, 774, 774,
			774, 775, 775, 776, 776, 777, 777, 778, 778, 778, 779, 779, 780, 780, 781, 781, 781, 782, 782, 783, 783,
			784, 784, 784, 785, 785, 786, 786, 787, 787, 788, 788, 788, 789, 789, 790, 790, 791, 791, 791, 792, 792,
			793, 793, 794, 794, 794, 795, 795, 796, 796, 797, 797, 798, 798, 798, 799, 799, 800, 800, 801, 801, 801,
			802, 802, 803, 803, 804, 804, 804, 805, 805, 806, 806, 807, 807, 808, 808, 808, 809, 809, 810, 810, 811,
			811, 811, 812, 812, 813, 813, 814, 814, 815, 815, 815, 816, 816, 817, 817, 818 };

	// sum of data criteria for pm = 0.75, retrieved from
	// https://github.com/Benson-Genomics-Lab/TRF/blob/master/src/tr30dat.c
	int sumdata75[] = { 0, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 6, 6,
			7, 7, 8, 8, 8, 9, 9, 10, 10, 11, 11, 11, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 10, 11, 11, 11, 12, 12, 12, 13,
			13, 13, 14, 14, 14, 15, 15, 16, 16, 16, 17, 17, 17, 18, 18, 18, 19, 19, 20, 20, 20, 21, 21, 21, 22, 22, 23,
			23, 23, 24, 24, 24, 25, 25, 26, 26, 26, 27, 27, 27, 28, 28, 29, 29, 29, 30, 30, 30, 31, 31, 32, 32, 32, 33,
			33, 34, 34, 34, 35, 35, 35, 36, 36, 37, 37, 37, 38, 38, 39, 39, 39, 40, 40, 41, 41, 41, 42, 42, 42, 43, 43,
			44, 44, 44, 45, 45, 46, 46, 46, 47, 47, 48, 48, 48, 49, 49, 50, 24, 24, 24, 24, 25, 25, 25, 25, 26, 26, 26,
			26, 27, 27, 27, 27, 28, 28, 28, 28, 29, 29, 29, 29, 30, 30, 30, 30, 31, 31, 31, 31, 32, 32, 32, 32, 33, 33,
			33, 34, 34, 34, 34, 35, 35, 35, 35, 36, 36, 36, 36, 37, 37, 37, 37, 38, 38, 38, 38, 39, 39, 39, 39, 40, 40,
			40, 40, 41, 41, 41, 41, 42, 42, 42, 43, 43, 43, 43, 44, 44, 44, 44, 45, 45, 45, 45, 46, 46, 46, 46, 47, 47,
			47, 48, 48, 48, 48, 49, 49, 49, 49, 50, 50, 50, 50, 51, 51, 51, 51, 52, 52, 52, 53, 53, 53, 53, 54, 54, 54,
			54, 55, 55, 55, 55, 56, 56, 56, 57, 57, 57, 57, 58, 58, 58, 58, 59, 59, 59, 60, 60, 60, 60, 61, 61, 61, 61,
			62, 62, 62, 62, 63, 63, 63, 64, 64, 64, 64, 65, 65, 65, 65, 66, 66, 66, 67, 67, 67, 67, 68, 68, 68, 68, 69,
			69, 69, 70, 70, 70, 70, 71, 71, 71, 71, 72, 72, 72, 73, 73, 73, 73, 74, 74, 74, 74, 75, 75, 75, 76, 76, 76,
			76, 77, 77, 77, 78, 78, 78, 78, 79, 79, 79, 79, 80, 80, 80, 81, 81, 81, 81, 82, 82, 82, 82, 83, 83, 83, 84,
			84, 84, 84, 85, 85, 85, 86, 86, 86, 86, 87, 87, 87, 87, 88, 88, 88, 89, 89, 89, 89, 90, 90, 90, 91, 91, 91,
			91, 92, 92, 92, 92, 93, 93, 93, 94, 94, 94, 94, 95, 95, 95, 96, 96, 96, 96, 97, 97, 97, 98, 98, 98, 98, 99,
			99, 99, 99, 100, 100, 100, 101, 101, 101, 101, 102, 102, 102, 103, 103, 103, 103, 104, 104, 104, 105, 105,
			105, 105, 106, 106, 106, 107, 107, 107, 107, 108, 108, 108, 109, 109, 109, 109, 110, 110, 110, 110, 111,
			111, 111, 112, 112, 112, 112, 113, 113, 113, 114, 114, 114, 114, 115, 115, 115, 116, 116, 116, 116, 117,
			117, 117, 118, 118, 118, 118, 119, 119, 119, 120, 120, 120, 120, 121, 121, 121, 122, 122, 122, 122, 123,
			123, 123, 124, 124, 124, 124, 125, 125, 125, 126, 126, 126, 126, 127, 127, 127, 128, 128, 128, 128, 129,
			129, 129, 130, 130, 130, 130, 131, 131, 131, 132, 132, 132, 132, 133, 133, 133, 134, 134, 134, 134, 135,
			135, 135, 136, 136, 136, 136, 137, 137, 137, 138, 138, 138, 138, 139, 139, 139, 140, 140, 140, 140, 141,
			141, 141, 142, 142, 142, 143, 143, 143, 143, 144, 144, 144, 145, 145, 145, 145, 146, 146, 146, 147, 147,
			147, 147, 148, 148, 148, 149, 149, 149, 149, 150, 150, 150, 151, 151, 151, 151, 152, 152, 152, 153, 153,
			153, 153, 154, 154, 154, 155, 155, 155, 156, 156, 156, 156, 157, 157, 157, 158, 158, 158, 158, 159, 159,
			159, 160, 160, 160, 160, 161, 161, 161, 162, 162, 162, 162, 163, 163, 163, 164, 164, 164, 165, 165, 165,
			165, 166, 166, 166, 167, 167, 167, 167, 168, 168, 168, 169, 169, 169, 169, 170, 170, 170, 171, 171, 171,
			172, 172, 172, 172, 173, 173, 173, 174, 174, 174, 174, 175, 175, 175, 176, 176, 176, 176, 177, 177, 177,
			178, 178, 178, 179, 179, 179, 179, 180, 180, 180, 181, 181, 181, 181, 182, 182, 182, 183, 183, 183, 183,
			184, 184, 184, 185, 185, 185, 186, 186, 186, 186, 187, 187, 187, 188, 188, 188, 188, 189, 189, 189, 190,
			190, 190, 191, 191, 191, 191, 192, 192, 192, 193, 193, 193, 193, 194, 194, 194, 195, 195, 195, 196, 196,
			196, 196, 197, 197, 197, 198, 198, 198, 198, 199, 199, 199, 200, 200, 200, 201, 201, 201, 201, 202, 202,
			202, 203, 203, 203, 203, 204, 204, 204, 205, 205, 205, 206, 206, 206, 206, 207, 207, 207, 208, 208, 208,
			208, 209, 209, 209, 210, 210, 210, 211, 211, 211, 211, 212, 212, 212, 213, 213, 213, 213, 214, 214, 214,
			215, 215, 215, 216, 216, 216, 216, 217, 217, 217, 218, 218, 218, 219, 219, 219, 219, 220, 220, 220, 221,
			221, 221, 221, 222, 222, 222, 223, 223, 223, 224, 224, 224, 224, 225, 225, 225, 226, 226, 226, 227, 227,
			227, 227, 228, 228, 228, 229, 229, 229, 229, 230, 230, 230, 231, 231, 231, 232, 232, 232, 232, 233, 233,
			233, 234, 234, 234, 235, 235, 235, 235, 236, 236, 236, 237, 237, 237, 237, 238, 238, 238, 239, 239, 239,
			240, 240, 240, 240, 241, 241, 241, 242, 242, 242, 243, 243, 243, 243, 244, 244, 244, 245, 245, 245, 246,
			246, 246, 246, 247, 247, 247, 248, 248, 248, 248, 249, 249, 249, 250, 250, 250, 251, 251, 251, 251, 252,
			252, 252, 253, 253, 253, 254, 254, 254, 254, 255, 255, 255, 256, 256, 256, 257, 257, 257, 257, 258, 258,
			258, 259, 259, 259, 260, 260, 260, 260, 261, 261, 261, 262, 262, 262, 263, 263, 263, 263, 264, 264, 264,
			265, 265, 265, 265, 266, 266, 266, 267, 267, 267, 268, 268, 268, 268, 269, 269, 269, 270, 270, 270, 271,
			271, 271, 271, 272, 272, 272, 273, 273, 273, 274, 274, 274, 274, 275, 275, 275, 276, 276, 276, 277, 277,
			277, 277, 278, 278, 278, 279, 279, 279, 280, 280, 280, 280, 281, 281, 281, 282, 282, 282, 283, 283, 283,
			283, 284, 284, 284, 285, 285, 285, 286, 286, 286, 286, 287, 287, 287, 288, 288, 288, 289, 289, 289, 289,
			290, 290, 290, 291, 291, 291, 292, 292, 292, 292, 293, 293, 293, 294, 294, 294, 295, 295, 295, 295, 296,
			296, 296, 297, 297, 297, 298, 298, 298, 298, 299, 299, 299, 300, 300, 300, 301, 301, 301, 301, 302, 302,
			302, 303, 303, 303, 304, 304, 304, 304, 305, 305, 305, 306, 306, 306, 307, 307, 307, 307, 308, 308, 308,
			309, 309, 309, 310, 310, 310, 311, 311, 311, 311, 312, 312, 312, 313, 313, 313, 314, 314, 314, 314, 315,
			315, 315, 316, 316, 316, 317, 317, 317, 317, 318, 318, 318, 319, 319, 319, 320, 320, 320, 320, 321, 321,
			321, 322, 322, 322, 323, 323, 323, 323, 324, 324, 324, 325, 325, 325, 326, 326, 326, 326, 327, 327, 327,
			328, 328, 328, 329, 329, 329, 330, 330, 330, 330, 331, 331, 331, 332, 332, 332, 333, 333, 333, 333, 334,
			334, 334, 335, 335, 335, 336, 336, 336, 336, 337, 337, 337, 338, 338, 338, 339, 339, 339, 339, 340, 340,
			340, 341, 341, 341, 342, 342, 342, 343, 343, 343, 343, 344, 344, 344, 345, 345, 345, 346, 346, 346, 346,
			347, 347, 347, 348, 348, 348, 349, 349, 349, 349, 350, 350, 350, 351, 351, 351, 352, 352, 352, 353, 353,
			353, 353, 354, 354, 354, 355, 355, 355, 356, 356, 356, 356, 357, 357, 357, 358, 358, 358, 359, 359, 359,
			359, 360, 360, 360, 361, 361, 361, 362, 362, 362, 363, 363, 363, 363, 364, 364, 364, 365, 365, 365, 366,
			366, 366, 366, 367, 367, 367, 368, 368, 368, 369, 369, 369, 369, 370, 370, 370, 371, 371, 371, 372, 372,
			372, 373, 373, 373, 373, 374, 374, 374, 375, 375, 375, 376, 376, 376, 376, 377, 377, 377, 378, 378, 378,
			379, 379, 379, 380, 380, 380, 380, 381, 381, 381, 382, 382, 382, 383, 383, 383, 383, 384, 384, 384, 385,
			385, 385, 386, 386, 386, 387, 387, 387, 387, 388, 388, 388, 389, 389, 389, 390, 390, 390, 390, 391, 391,
			391, 392, 392, 392, 393, 393, 393, 394, 394, 394, 394, 395, 395, 395, 396, 396, 396, 397, 397, 397, 397,
			398, 398, 398, 399, 399, 399, 400, 400, 400, 401, 401, 401, 401, 402, 402, 402, 403, 403, 403, 404, 404,
			404, 404, 405, 405, 405, 406, 406, 406, 407, 407, 407, 408, 408, 408, 408, 409, 409, 409, 410, 410, 410,
			411, 411, 411, 411, 412, 412, 412, 413, 413, 413, 414, 414, 414, 415, 415, 415, 415, 416, 416, 416, 417,
			417, 417, 418, 418, 418, 419, 419, 419, 419, 420, 420, 420, 421, 421, 421, 422, 422, 422, 422, 423, 423,
			423, 424, 424, 424, 425, 425, 425, 426, 426, 426, 426, 427, 427, 427, 428, 428, 428, 429, 429, 429, 430,
			430, 430, 430, 431, 431, 431, 432, 432, 432, 433, 433, 433, 433, 434, 434, 434, 435, 435, 435, 436, 436,
			436, 437, 437, 437, 437, 438, 438, 438, 439, 439, 439, 440, 440, 440, 441, 441, 441, 441, 442, 442, 442,
			443, 443, 443, 444, 444, 444, 444, 445, 445, 445, 446, 446, 446, 447, 447, 447, 448, 448, 448, 448, 449,
			449, 449, 450, 450, 450, 451, 451, 451, 452, 452, 452, 452, 453, 453, 453, 454, 454, 454, 455, 455, 455,
			456, 456, 456, 456, 457, 457, 457, 458, 458, 458, 459, 459, 459, 459, 460, 460, 460, 461, 461, 461, 462,
			462, 462, 463, 463, 463, 463, 464, 464, 464, 465, 465, 465, 466, 466, 466, 467, 467, 467, 467, 468, 468,
			468, 469, 469, 469, 470, 470, 470, 471, 471, 471, 471, 472, 472, 472, 473, 473, 473, 474, 474, 474, 475,
			475, 475, 475, 476, 476, 476, 477, 477, 477, 478, 478, 478, 479, 479, 479, 479, 480, 480, 480, 481, 481,
			481, 482, 482, 482, 483, 483, 483, 483, 484, 484, 484, 485, 485, 485, 486, 486, 486, 486, 487, 487, 487,
			488, 488, 488, 489, 489, 489, 490, 490, 490, 490, 491, 491, 491, 492, 492, 492, 493, 493, 493, 494, 494,
			494, 494, 495, 495, 495, 496, 496, 496, 497, 497, 497, 498, 498, 498, 498, 499, 499, 499, 500, 500, 500,
			501, 501, 501, 502, 502, 502, 502, 503, 503, 503, 504, 504, 504, 505, 505, 505, 506, 506, 506, 506, 507,
			507, 507, 508, 508, 508, 509, 509, 509, 510, 510, 510, 510, 511, 511, 511, 512, 512, 512, 513, 513, 513,
			514, 514, 514, 514, 515, 515, 515, 516, 516, 516, 517, 517, 517, 518, 518, 518, 518, 519, 519, 519, 520,
			520, 520, 521, 521, 521, 522, 522, 522, 522, 523, 523, 523, 524, 524, 524, 525, 525, 525, 526, 526, 526,
			526, 527, 527, 527, 528, 528, 528, 529, 529, 529, 530, 530, 530, 530, 531, 531, 531, 532, 532, 532, 533,
			533, 533, 534, 534, 534, 534, 535, 535, 535, 536, 536, 536, 537, 537, 537, 538, 538, 538, 539, 539, 539,
			539, 540, 540, 540, 541, 541, 541, 542, 542, 542, 543, 543, 543, 543, 544, 544, 544, 545, 545, 545, 546,
			546, 546, 547, 547, 547, 547, 548, 548, 548, 549, 549, 549, 550, 550, 550, 551, 551, 551, 551, 552, 552,
			552, 553, 553, 553, 554, 554, 554, 555, 555, 555, 555, 556, 556, 556, 557, 557, 557, 558, 558, 558, 559,
			559, 559, 559, 560, 560, 560, 561, 561, 561, 562, 562, 562, 563, 563, 563, 564, 564, 564, 564, 565, 565,
			565, 566, 566, 566, 567 };

	List<Integer> primeNumbers;

	/**
	 * Initializes this instance of a candidate selector
	 * 
	 * @param sequence
	 *            sequence to be analyzed
	 * @param sequenceName
	 *            name of the sequence
	 * @param matchProbability
	 *            can be 0.8 or 0.75
	 * @param indelProbability
	 *            typically 0.1
	 * @param maxSize
	 *            maximum size of the pattern
	 * @param minimumAlignmentScore
	 *            minimum score for a pattern to be reported
	 */
	public TRFCandidateSelector(String sequence, String sequenceName, double matchProbability, double indelProbability,
			int maxSize, int minimumAlignmentScore, int matchScore, int missScore) {

		// initializes empty atributes
		this.indelProbability = indelProbability;
		this.matchProbability = matchProbability;
		this.sequenceName = sequenceName;
		this.sequence = sequence;
		this.matchScore = matchScore;
		this.missScore = missScore;
		this.historyList = new HashMap<String, ArrayList<Integer>>();
		this.distanceList = new HashMap<Integer, ArrayList<Integer>>();
		this.allCandidates = new ArrayList<TandemRepeat>();
		this.minimumAlignmentScore = minimumAlignmentScore;
		int[] tupleSizes = { 5 };
		this.maxPatternSize = maxSize;
		HashMap<Integer, TreeMap<Integer, ArrayList<TandemRepeat>>> initialC = null;
		ArrayList<TandemRepeat> candidate = null;

		int mid = maxSize / 2;
		this.primeNumbers = primeNumbersBruteForce(mid);
		long time = System.currentTimeMillis();

		// in case of multiple tuple sizes, candidates are searched for each
		// tuple size
		for (int i : tupleSizes) {
			initialC = getCandidates(i);
			candidate = refineCandidates(initialC);
			candidate.sort((o1, o2) -> Integer.compare(o1.getFirst(), o2.getFirst()));
			removeOverlaps(candidate);
			this.allCandidates.addAll(candidate);
			initialC = null;
			candidate = null;

		}
		// after finding the candidates, overlap tandem repeats are removed
		this.allCandidates.sort((o1, o2) -> Integer.compare(o1.getFirst(), o2.getFirst()));
		removeOverlaps(this.allCandidates);
		time = System.currentTimeMillis() - time;
		System.out.println("Time searching candidates: " + time + "ms");

	}

	/**
	 * Method that finds candidates using TRF definition mainly
	 * 
	 * @param kTupleSize
	 * @return hashmap containing as keys distance of the TR, values another
	 *         treemap containing index where the match was reported and list of
	 *         TR at that index
	 */
	public HashMap<Integer, TreeMap<Integer, ArrayList<TandemRepeat>>> getCandidates(int kTupleSize) {
		// Answer structure
		HashMap<Integer, TreeMap<Integer, ArrayList<TandemRepeat>>> candidates = new HashMap<Integer, TreeMap<Integer, ArrayList<TandemRepeat>>>();
		ArrayList<Integer> occurrences = null;
		ArrayList<Integer> matchesAtD = null;
		ArrayList<TandemRepeat> tra = null;
		TreeMap<Integer, ArrayList<TandemRepeat>> indexes = null;
		// loop over the sequence to find candidates
		for (int i = kTupleSize - 1; i < sequence.length(); i++) {
			// i final index of probe p, inclusive
			String probe = sequence.substring(i - kTupleSize + 1, i + 1);
			// check previous occurrences of this probe
			occurrences = historyList.getOrDefault(probe, new ArrayList<Integer>());
			// now we scan for occurrences of probe
			if (occurrences.size() > 0) { // contains occurrences of this probe
				for (int u = occurrences.size() - 1; u >= 0; u--) {
					int j = occurrences.get(u);
					int d = i - j;
					// If d is greater than max size we will not calculate this
					// tandem repeat
					if (d > maxPatternSize) {
						break;
					} else {
						int totalMatches = 0;
						int leftMostIndex = Integer.MAX_VALUE;
						// Max delta distance to the left and right to check for
						// indels
						double maxDeltaDistance = 2.3 * Math.sqrt(indelProbability * d);
						int maxDelta = (int) Math.floor(maxDeltaDistance);
						// We get matches at this distance d
						matchesAtD = distanceList.getOrDefault(d, new ArrayList<Integer>());
						// We add this match index to the list
						matchesAtD.add(i);
						int maxRange = i - d + 1;
						// matches that happened before i-d+1 are removed, d
						// could be taken as a sliding window
						matchesAtD.removeIf(s -> s < maxRange);
						// We update distances list with new matches
						distanceList.put(d, matchesAtD);
						// Now we search for the minimum element in matches
						// array
						Integer min = matchesAtD.stream().mapToInt(v -> v).min().orElse(Integer.MAX_VALUE) - kTupleSize
								+ 1; // moderate
						if (min < leftMostIndex) {
							leftMostIndex = min;
						}
						// we calculate sum of heads for this probe
						totalMatches += matchesAtD.size() * kTupleSize;
						// update nearby distances list
						for (int k = 1; k < maxDelta + 1; k++) {
							int dleft = d - k;
							int dright = d + k;
							// if a match was found at left
							if (dleft > 0) {
								// perform same operations for d
								if (distanceList.containsKey(dleft)) {
									matchesAtD = distanceList.get(dleft);
									matchesAtD.removeIf(s -> s < maxRange);
									distanceList.put(dleft, matchesAtD);
									totalMatches += matchesAtD.size() * kTupleSize;
									min = matchesAtD.stream().mapToInt(v -> v).min().orElse(Integer.MAX_VALUE)
											- kTupleSize + 1;
									if (min < leftMostIndex) {
										leftMostIndex = min;
									}
								}
							}
							// if a match was found at right
							if (dright > 0) {
								// perform same operations for d
								if (distanceList.containsKey(dright)) {
									matchesAtD = distanceList.get(dright);
									matchesAtD.removeIf(s -> s < maxRange);
									distanceList.put(dright, matchesAtD);
									totalMatches += matchesAtD.size() * kTupleSize;
									min = matchesAtD.stream().mapToInt(v -> v).min().orElse(Integer.MAX_VALUE)
											- kTupleSize + 1;
									if (min < leftMostIndex) {
										leftMostIndex = min;
									}
								}
							}
						}
						// now we compute apparent size to test criteria
						int apparentSize = i - leftMostIndex + 1;
						// limit apparent size
						if (apparentSize > d) {
							apparentSize = d;
						}
						// compute and test criteria
						int apTh = computeApparentSizeThreshold(d);
						int suOh = computeSumOfHeadsThreshold(d, kTupleSize);
						if (totalMatches >= suOh && apparentSize >= apTh) {
							// if TR matches criteria it is created
							TandemRepeat trc = new TandemRepeat(j + 1, i, apparentSize, totalMatches);
							indexes = candidates.getOrDefault(apparentSize,
									new TreeMap<Integer, ArrayList<TandemRepeat>>());
							// new candidate is added in index i
							tra = indexes.getOrDefault(i, new ArrayList<TandemRepeat>());
							tra.add(trc);
							indexes.put(i, tra);
							// tr with apparent size is added to the answer
							// hashmap
							candidates.put(apparentSize, indexes);
						}
					}
				}
			}
			// finally we update historyList list
			occurrences.add(i);
			historyList.put(probe, occurrences);
		}

		return candidates;
	}

	/**
	 * This method remove redundant candidates
	 * 
	 * @param candidates
	 * @return
	 */
	public ArrayList<TandemRepeat> refineCandidates(
			HashMap<Integer, TreeMap<Integer, ArrayList<TandemRepeat>>> candidates) {
		ArrayList<TandemRepeat> answer = new ArrayList<TandemRepeat>();

		// refine candidates of same apparent size
		for (Integer apSize : candidates.keySet()) {
			boolean notEnough = true;
			while (notEnough) {
				// matches with apsize d at index i
				TreeMap<Integer, ArrayList<TandemRepeat>> indexes = candidates.get(apSize);

				double maxDeltaDistance = 2.3 * Math.sqrt(indelProbability * apSize);
				// maximum size at which there could be a match
				int maxDelta = (int) Math.floor(maxDeltaDistance);
				boolean search = true;
				int original = indexes.lastKey(); // we know original index

				int largest = indexes.lastKey(); // we know original index
				int numCopies = 1;
				ArrayList<TandemRepeat> miniCandi = indexes.get(original);
				TandemRepeat finalC = miniCandi.get(0);
				int lowerborder = largest - apSize - maxDelta;
				int upperborder = largest - apSize;
				indexes.remove(largest);
				int beginning = finalC.getPatternBeginning() - finalC.getDistance(); // default
																						// beginning
				while (search) {
					// i dint find more redundant TR at this distance
					if (indexes.keySet().isEmpty()) {
						search = false;
						notEnough = false; // is enough
					} else {
						largest = indexes.lastKey(); // next one
						if (largest > upperborder) { // this Tr is redundant if
														// it is above the
														// border
							indexes.remove(largest);
						} else if (largest >= lowerborder && largest <= upperborder) { // count,
																						// found
																						// another
																						// occurence
							// compute beginning
							beginning = indexes.get(largest).get(0).getPatternBeginning() - apSize;
							// update num copies
							numCopies++;
							// update border
							lowerborder = largest - apSize - maxDelta;
							upperborder = largest - apSize;

						} else { // there is another TR at lower index
							candidates.put(apSize, indexes);
							break;
						}

					}

				}
				// for definition if a TR pass criteria num copies is at least
				// of 2
				numCopies++;
				// update TR and add it to answer
				finalC.setNumCopies(numCopies);
				finalC.setFirst(beginning);
				answer.add(finalC);

			}

		}

		return answer;
	}

	/**
	 * Method that removes tandem repeats of different size and num copies
	 * overlapping the same section
	 * 
	 * @param candidates
	 */
	public void removeOverlaps(ArrayList<TandemRepeat> candidates) {
		if (candidates.size() == 0)
			return;
		TandemRepeat bestC = candidates.get(0);
		int indexBest = 0;
		// minimum index of section
		int left = bestC.getFirst();
		// maximum index of section
		int right = bestC.getLast();
		int size = bestC.getTotalSize();
		for (int i = 1; i < candidates.size(); i++) {
			TandemRepeat trc = candidates.get(i);
			int begin = trc.getFirst();
			int end = trc.getLast();
			// if the 2 TR are overlapped we update ref indexes
			if (isOverlapped(left, right, begin, end)) {
				if (begin < left) {
					left = begin;
				}
				if (end > right) {
					right = end;
				}
				if (trc.getTotalSize() >= size) { // best, covering more
													// distance
					// update best candidate
					candidates.remove(indexBest);
					size = trc.getTotalSize();
					i--;
					indexBest = i;
				} else {
					candidates.remove(i);
					i--;
				}
			} else {
				// this is a different TR, we reset ref variables
				bestC = candidates.get(i);
				indexBest = i;
				left = bestC.getFirst();
				right = bestC.getLast();
				size = bestC.getTotalSize();
			}

		}
	}

	/**
	 * This method allign candidates with ideal sequence and refine num copies
	 * 
	 * @param candidates
	 * @throws IOException
	 */
	public void allignCandidates(ArrayList<TandemRepeat> candidates) throws IOException {

		SimpleEditDistanceMeasure a = new SimpleEditDistanceMeasure();
		PrintWriter pw = new PrintWriter(new FileWriter("data\\outTR.txt", true));
		String line = "";
		for (TandemRepeat trCandidate : candidates) {
			// we compute pattern
			String pattern = sequence.substring(trCandidate.getPatternBeginning(), trCandidate.getLast() + 1);
			double times = trCandidate.getNumCopies();

			String realseq = sequence.substring(trCandidate.getFirst(), trCandidate.getLast() + 1);

			String[] result = cutPattern(pattern, times);
			pattern = result[0];
			times = Double.parseDouble(result[1]);

			// extend sequence to get times
			int beg = trCandidate.getFirst();
			int patternSize = pattern.length();
			List<CharSequence> temp = null;
			boolean keepGoing = true;

			// extend backwards N times in case algorithm didn't find these
			// copies
			beg -= patternSize;
			while (keepGoing) {

				if (beg >= 0) {
					temp = a.pairwiseAlignment(pattern, sequence.substring(beg, beg + patternSize));
					String s1 = temp.get(0).toString();
					String s2 = temp.get(1).toString();
					int ham = hammingDist(s1, s2);

					double limit = 0.4 * s1.length();

					if (ham <= limit) {
						times++;
						beg -= patternSize;
					} else {
						keepGoing = false;
					}
				} else {
					keepGoing = false;
				}

			}
			beg += patternSize;
			int beginning = trCandidate.getFirst();
			int reallength = realseq.length();

			// if we found more matching copies update info that is going to be
			// written
			if (times > trCandidate.getNumCopies()) {
				beginning = beg;
				reallength = trCandidate.getLast() - beg + 1;
				realseq = sequence.substring(beg, trCandidate.getLast() + 1);
			}

			// extend sequence forward
			int end = trCandidate.getLast() + 1;
			temp = null;
			double forwardtimes = times;
			keepGoing = true;

			// end += patternSize;
			while (keepGoing) {

				if (end + patternSize < sequence.length()) {
					temp = a.pairwiseAlignment(pattern, sequence.substring(end, end + patternSize));
					String s1 = temp.get(0).toString();
					String s2 = temp.get(1).toString();
					int ham = hammingDist(s1, s2);

					double limit = 0.4 * s1.length();

					if (ham <= limit) {
						forwardtimes++;
						end += patternSize;
					} else {
						keepGoing = false;
					}
				} else {
					keepGoing = false;
				}

			}
			// end -= patternSize;
			// if we found more matching copies forward we update info that is
			// going to be written
			if (forwardtimes > times) {
				times = forwardtimes;
				trCandidate.setLast(end);
				reallength = trCandidate.getLast() - beg + 1;
				realseq = sequence.substring(beg, trCandidate.getLast() + 1);
			}

			// calculate ideal seq by multiplying the pattern
			String toallign = "";
			for (int i = 0; i < times; i++) {
				toallign += pattern;
			}

			// allign final sequences
			List<CharSequence> r = a.pairwiseAlignment(toallign, realseq);

			// now we check if the TR is gonna be reported
			int alscore = alignmentScore(r.get(0).toString(), r.get(1).toString());
			boolean report = alscore >= minimumAlignmentScore;

			trCandidate.setSequenceName(sequenceName);

			// write line to file
			if (report && realseq.length() >= 25) {
				if ((times < 3.0 && pattern.length() >= 12) || times >= 3.0) {
					if ((times < 3.0 && alscore >= minimumAlignmentScore + 10)
							|| (times < 4.0 && times >= 3.0 && alscore >= minimumAlignmentScore + 2) || times >= 4.0) {
						line = sequenceName + " " + beginning + " " + trCandidate.getLast() + " "
								+ trCandidate.getDistance() + " " + times + " " + reallength + " " + pattern + " "
								+ realseq + " " + toallign + " " + r + " " + alscore;
						pw.println(line);
					}

				}
			}

		}
		pw.close();
	}

	public String[] cutPattern(String pattern, double times) {
		String[] ans = new String[2];

		String newpattern = pattern;
		double newtimes = times;

		ans[0] = pattern;
		ans[1] = times + "";
		int patLength = newpattern.length();
		for (Integer n : primeNumbers) {
			// if a pattern can be split in two, we split it until necessary
			// it can be split by this number
			if (n == 2) {
				int mid = newpattern.length() / 2; // get the middle of the
													// String
				boolean keepCutting = true;
				while (keepCutting) {
					String[] parts = { newpattern.substring(0, mid), newpattern.substring(mid) };
					if (parts[0].equals(parts[1])) {
						newpattern = parts[0];
						newtimes = times * 2;
						mid = newpattern.length() / 2;
					} else {
						keepCutting = false;
					}
				}
				if (!newpattern.equals(pattern)) {
					ans[0] = newpattern;
					ans[1] = newtimes + "";
					return ans;
				}
			} else { // any other prime number
				if (patLength % n == 0) {
					boolean couldCut = true;
					ArrayList<String> subpatterns = new ArrayList<String>();
					int limit = 0;
					for (int i = 0; i < patLength; i += n) {
						limit = i + n;
						String subpattern = pattern.substring(i, limit);
						subpatterns.add(subpattern);
					}
					String firstp = subpatterns.get(0);
					for (int i = 1; i < subpatterns.size(); i++) {
						couldCut &= firstp.equals(subpatterns.get(i));
						firstp = subpatterns.get(i);
						if (!couldCut) {
							break;
						}
					}
					if (couldCut) {
						newpattern = firstp;
						newtimes = newtimes * n;
						ans[0] = newpattern;
						ans[1] = newtimes + "";
						return ans;
					}
				}
			}

		}

		return ans;
	}

	public static List<Integer> primeNumbersBruteForce(int n) {
		List<Integer> primeNumbers = new LinkedList<>();
		for (int i = 2; i <= n; i++) {
			if (isPrimeBruteForce(i)) {
				primeNumbers.add(i);
			}
		}
		return primeNumbers;
	}

	public static boolean isPrimeBruteForce(int number) {
		for (int i = 2; i < number; i++) {
			if (number % i == 0) {
				return false;
			}
		}
		return true;
	}

	/*
	 * Computes hamming distance from two strings
	 * 
	 * @param str1
	 * 
	 * @param str2
	 * 
	 * @return
	 */
	public int hammingDist(String str1, String str2) {
		int i = 0, count = 0;
		while (i < str1.length()) {
			if (str1.charAt(i) != str2.charAt(i))
				count++;
			i++;
		}
		return count;
	}

	/**
	 * Checks if two alligned strings have at least N same characters
	 * 
	 * @param str1
	 * @param str2
	 * @param minimum
	 * @return
	 */
	public int alignmentScore(String str1, String str2) {
		int i = 0, count = 0;

		while (i < str1.length()) {
			if (str1.charAt(i) == '-' || str2.charAt(i) == '-') {
				count -= missScore;
			} else if (str1.charAt(i) != str2.charAt(i)) {
				count -= missScore;
			} else if (str1.charAt(i) == str2.charAt(i))
				count += matchScore;
			i++;
		}
		return count;
	}

	/**
	 * Checks if sequence from i1 to i2 overlaps with sequence from j1 to j2
	 * 
	 * @param i1
	 * @param i2
	 * @param j1
	 * @param j2
	 * @return
	 */
	public boolean isOverlapped(int i1, int i2, int j1, int j2) {
		boolean ans = false;
		ans |= (i1 == j1 && i2 == j2);
		ans |= (j1 <= i1 && (j2 >= i1 && j2 <= i2));
		ans |= ((j1 >= i1 && j1 <= i2) && j2 >= i2);
		ans |= (j1 <= i1 && j2 >= i2);
		ans |= (j1 >= i1 && j2 <= i2);
		return ans;
	}

	/**
	 * Returns sum of heads criteria for distance d and tuple size k
	 * 
	 * @param d
	 * @param k
	 * @return
	 */
	public int computeSumOfHeadsThreshold(int d, int k) {
		// return sumdata75[d - 1];
		if (matchProbability == 0.8) {
			return sumdata80[d - 1];
		}
		return sumdata75[d - 1];
	}

	/**
	 * Returns sum of heads criteria for distance d
	 * 
	 * @param d
	 * @param k
	 * @return
	 */
	public int computeApparentSizeThreshold(int d) {
		if (d > 10) {
			return (int) Math.floor(d * 0.9);
		} else if (d > 5) {
			return (int) Math.floor(d * 0.85);
		} else {
			return (int) Math.floor(d * 0.8);
		}

	}

	/**
	 * @return candidates found
	 */
	public ArrayList<TandemRepeat> getAllCandidates() {
		return allCandidates;
	}

	/**
	 * @return the sequence
	 */
	public String getSequence() {
		return sequence;
	}

	/**
	 * @return the historyList
	 */
	public HashMap<String, ArrayList<Integer>> getHistoryList() {
		return historyList;
	}

	/**
	 * @return the matchProbability
	 */
	public double getMatchProbability() {
		return matchProbability;
	}

	/**
	 * @return the indelProbability
	 */
	public double getIndelProbability() {
		return indelProbability;
	}

}
