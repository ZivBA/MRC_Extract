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

public class Willison64PossibilitiesYaoLike extends MeshiProgram implements Residues,AtomTypes {

	public static void main(String[] args) throws Exception {
//		String[] runs = {"ZQHBEDAG ZQHBEDAG 12_6784_7156",
//				"ZQHEGDAB ZQHEGDAB 12_6668_7048",
//				"ZGHEBDAQ ZGHEBDAQ 12_6618_6980",
//				"QZHEBDAG QZHEBDAG 12_6578_6932",
//				"ZQHDBGAE ZQHDBGAE 10_6026_6588",
//				"ZBQEHDAG ZBQEHDAG 10_5956_6514",
//				"HQBEZDAG HQBEZDAG 10_5912_6478",
//				"ZHQDBEAG ZHQDBEAG 08_5362_6108",
//				"ZDHEGQAB ZDHEGQAB 08_5162_5900",
//				"ZHQEBDAG GZHQEBDA 07_4841_5655",
//				"ZDGEQHAB EQHABZDG 07_4819_5643",
//				"ZGHDBQAE BQAEZGHD 06_4500_5421",
//				"GHDEBQAZ GHDEBQAZ 06_4400_5312",
//				"BDHGQEAZ QEAZBDHG 05_4084_5087",
//				"HQZEDBAG AGHQZEDB 05_4054_5055",
//				"ZEDQHGAB ZEDQHGAB 04_3748_4846",
//				"ZBHQDGAE DGAEZBHQ 04_3716_4813",
//				"EDBGHZAQ ZAQEDBGH 03_3396_4587",
//				"DQZEGHAB BDQZEGHA 03_3348_4535",
//				"QBDHEZAG HEZAGQBD 02_2972_4252",
//				"QEBZHGAD EBZHGADQ 01_2579_3954"};

		String[] runs = {"HZQBDEAG ZQBDEAGH 04_3888_4969",
				"QGHEDBAZ QGHEDBAZ 06_4398_5292",
				"HGZQEBAD QEBADHGZ 01_2615_3984",
				"GZEBDHAQ AQGZEBDH 02_3022_4294",
				"ZHEBDGAQ GAQZHEBD 02_3079_4356",
				"ZQBEDHAG ZQBEDHAG 10_6020_6564",
				"GEBDHZAQ HZAQGEBD 01_2713_4083",
				"ZQEDHBAG GZQEDHBA 05_4211_5181",
				"BZQHGDAE DAEBZQHG 03_3350_4492",
				"ZQGBDHAE BDHAEZQG 05_4152_5162",
				"QZGEBDAH GEBDAHQZ 04_3760_4881",
				"DHZQGEAB EABDHZQG 02_3023_4296",
				"GHQEBDAZ HQEBDAZG 06_4445_5313",
				"EDHBZQAG ZQAGEDHB 06_4522_5416",
				"EBDHQGAZ DHQGAZEB 01_2610_3949",
				"DBHZGQAE DBHZGQAE 04_3630_4748",
				"QEBHGDAZ AZQEBHGD 04_3689_4777",
				"HBGZQDAE AEHBGZQD 03_3268_4503",
				"BDZEQHAG EQHAGBDZ 05_4040_5024",
				"QHZGBDAE GBDAEQHZ 03_3349_4509",
				"QZHGEBAD EBADQZHG 03_3380_4544",
				"GQEDHBAZ GQEDHBAZ 04_3850_4894",
				"EBQGZHAD EBQGZHAD 02_2910_4226",
				"EZBQDHAG QDHAGEZB 03_3356_4529",
				"ZQBDHGAE QBDHGAEZ 03_3451_4619",
				"GHZEBDAQ BDAQGHZE 04_3655_4761",
				"DBQGZHAE EDBQGZHA 01_2614_3989",
				"GZBDEHAQ GZBDEHAQ 02_3098_4340",
				"GEZHBDAQ QGEZHBDA 03_3325_4477",
				"GHEDBQAZ AZGHEDBQ 03_3403_4566",
				"QGZBEDAH AHQGZBED 02_2992_4268",
				"DHZGQEAB QEABDHZG 02_2998_4262",
				"ZQGDEBAH GDEBAHZQ 03_3483_4688",
				"GZHDBEAQ ZHDBEAQG 05_4183_5164",
				"EBDZQHAG GEBDZQHA 02_2991_4269",
				"GZEBDQAH GZEBDQAH 02_3018_4288",
				"BHEZDQAG ZDQAGBHE 03_3385_4545",
				"ZQHGDEAB QHGDEABZ 04_3871_4962",
				"HEDQGZAB HEDQGZAB 02_2890_4190",
				"EDQGHZAB BEDQGHZA 01_2585_3978",
				"DBHEQZAG QZAGDBHE 04_3703_4773",
				"ZQDHEBAG BAGZQDHE 05_4125_5126",
				"EZBDQGAH QGAHEZBD 01_2692_4059",
				"ZEBDQGAH AHZEBDQG 06_4452_5417",
				"EBDZHGAQ ZHGAQEBD 02_3044_4353",
				"QHZEBDAG DAGQHZEB 05_4056_5052",
				"EQZHGDAB GDABEQZH 03_3365_4539",
				"EDBHGQAZ HGQAZEDB 01_2628_3968",
				"GHQEDBAZ QEDBAZGH 02_3029_4293",
				"BEHDQGAZ DQGAZBEH 03_3382_4593",
				"BDEZGHAQ EZGHAQBD 01_2627_4004",
				"HGEZBDAQ GEZBDAQH 03_3296_4485",
				"ZGDEBHAQ HAQZGDEB 05_4041_5033",
				"BHDQGZAE QGZAEBHD 01_2615_3967",
				"HQBDGZAE BDGZAEHQ 02_3088_4376",
				"GBZQHEAD ZQHEADGB 06_4514_5444",
				"QBEDZHAG DZHAGQBE 03_3381_4536",
				"HGZBDQAE DQAEHGZB 03_3349_4546",
				"HDZQGBAE EHDZQGBA 01_2572_3970",
				"GDHBQZAE BQZAEGDH 03_3368_4563",
				"QGBHEDAZ DAZQGBHE 02_2987_4235",
				"ZBDHQGAE ZBDHQGAE 04_3704_4842",
				"BHQZEDAG QZEDAGBH 03_3440_4602",
				"GDZBQHAE ZBQHAEGD 02_3052_4340",
				"BGEDZQAH ZQAHBGED 04_3786_4896",
				"QGDBEHAZ EHAZQGDB 01_2617_3953",
				"HEDZBQAG BQAGHEDZ 04_3751_4846",
				"HBEDZQAG HBEDZQAG 04_3704_4794",
				"GQHDBEAZ HDBEAZGQ 05_4170_5168",
				"EGZQDBAH ZQDBAHEG 04_3786_4880",
				"EGQZDBAH BAHEGQZD 03_3273_4503",
				"ZBDQGHAE GHAEZBDQ 03_3348_4534",
				"ZBHDGQAE DGQAEZBH 03_3439_4628",
				"BHDQZEAG HDQZEAGB 02_2995_4263",
				"ZEGBDHAQ QZEGBDHA 04_3766_4858",
				"QEBDZHAG GQEBDZHA 03_3445_4588",
				"GQEZBDAH ZBDAHGQE 05_4110_5130",
				"QGDHEBAZ GDHEBAZQ 04_3682_4757",
				"HEDQZGAB DQZGABHE 02_2982_4286",
				"QGBHDEAZ BHDEAZQG 03_3398_4551",
				"ZDQGHEAB ZDQGHEAB 04_3860_4982",
				"ZHBDQEAG BDQEAGZH 04_3818_4952",
				"QEBGHZAD QEBGHZAD 02_3022_4336",
				"EZBHQGAD DEZBHQGA 01_2628_4013",
				"GHBDZEAQ EAQGHBDZ 01_2695_4049",
				"EBGDHZAQ EBGDHZAQ 02_3032_4338",
				"GBHDEZAQ QGBHDEZA 02_3050_4302",
				"EDZHGQAB ABEDZHGQ 01_2579_3946",
				"GQHZEBAD QHZEBADG 06_4435_5326",
				"EBHZDQAG HZDQAGEB 03_3308_4497",
				"GHDZQEAB ZQEABGHD 04_3773_4856",
				"ZGHQEDAB HQEDABZG 06_4520_5419",
				"GDBZQHAE HAEGDBZQ 01_2569_3955",
				"DZQHEBAG GDZQHEBA 02_3067_4314",
				"QGBDHEAZ DHEAZQGB 01_2718_4046",
				"GHBQZDAE BQZDAEGH 03_3457_4660",
				"DBGZHQAE HQAEDBGZ 03_3366_4564",
				"DEZQGBAH GBAHDEZQ 01_2606_3986",
				"QDEHGBAZ GBAZQDEH 02_2964_4201",
				"EQGBDZAH DZAHEQGB 02_3052_4315",
				"GEZHDBAQ BAQGEZHD 01_2605_3967",
				"BEHDZQAG QAGBEHDZ 03_3380_4577",
				"DZEHQGAB HQGABDZE 04_3686_4776",
				"DBZQHEAG QHEAGDBZ 03_3395_4533",
				"DEZBQGAH EZBQGAHD 01_2603_4003",
				"EZQDBHAG AGEZQDBH 04_3721_4808",
				"ZQEHDGAB ABZQEHDG 04_3726_4814",
				"EDHQZGAB QZGABEDH 03_3347_4572",
				"GQBEDZAH HGQBEDZA 04_3759_4843",
				"ZHGQBEAD QBEADZHG 04_3777_4868",
				"GQDBHEAZ AZGQDBHE 02_3070_4350",
				"EHQBDGAZ EHQBDGAZ 02_2994_4252",
				"EBGDQZAH ZAHEBGDQ 05_4072_5138",
				"ZHDBEQAG AGZHDBEQ 03_3389_4545",
				"BQGDEHAZ QGDEHAZB 03_3399_4584",
				"BGQZEHAD ADBGQZEH 01_2613_4006",
				"BDEQHZAG HZAGBDEQ 04_3710_4778",
				"BDGHZEAQ BDGHZEAQ 02_3012_4324",
				"HDEBZQAG GHDEBZQA 04_3688_4774",
				"GEHBZQAD HBZQADGE 03_3319_4541",
				"QHDEBZAG GQHDEBZA 06_4489_5378",
				"DQZHEGAB DQZHEGAB 04_3730_4800",
				"EZDQGBAH AHEZDQGB 01_2580_3958",
				"GHDBEZAQ DBEZAQGH 01_2663_4011",
				"QZGHEBAD ADQZGHEB 01_2608_3967",
				"GEBZHQAD ADGEBZHQ 03_3340_4571",
				"EHQBGZAD HQBGZADE 02_3005_4274",
				"GQZDEHAB DEHABGQZ 04_3730_4820",
				"BEQZHDAG EQZHDAGB 04_3710_4775",
				"GZDQHEAB BGZDQHEA 01_2661_4017",
				"BEGQZHAD ZHADBEGQ 03_3427_4660",
				"BHGZDQAE QAEBHGZD 01_2609_4019",
				"EHQBZDAG BZDAGEHQ 03_3390_4519",
				"BQGZHDAE BQGZHDAE 06_4486_5432",
				"GZDBQHAE ZDBQHAEG 03_3417_4593",
				"BDEQGHAZ BDEQGHAZ 02_2928_4206",
				"HEQGZBAD EQGZBADH 03_3324_4560",
				"BZHGDEAQ GDEAQBZH 02_2985_4269",
				"GBEZHQAD DGBEZHQA 02_2962_4248",
				"BQZEDHAG ZEDHAGBQ 05_4135_5145",
				"ZGEBHDAQ BHDAQZGE 03_3401_4556",
				"ZDBQEHAG DBQEHAGZ 04_3777_4866",
				"HGBZDEAQ AQHGBZDE 04_3726_4831",
				"GQEBHDAZ HDAZGQEB 03_3381_4518",
				"EGZBQHAD HADEGZBQ 02_2930_4216",
				"ZBDEGQAH BDEGQAHZ 03_3345_4553",
				"ZEQBDGAH QBDGAHZE 02_3026_4352",
				"HBEDQGAZ BEDQGAZH 01_2569_3979",
				"BDQEZGAH DQEZGAHB 03_3339_4524",
				"DGZHBQAE QAEDGZHB 02_2978_4228",
				"QBHDGZAE BHDGZAEQ 02_2983_4245",
				"HQGEZBAD HQGEZBAD 06_4418_5396",
				"QEGBZDAH DAHQEGBZ 03_3371_4597",
				"QHGZBEAD HGZBEADQ 02_2999_4278",
				"HBGQZDAE QZDAEHBG 03_3360_4541",
				"HQZDEGAB BHQZDEGA 02_3052_4323",
				"BZQEHGAD HGADBZQE 03_3351_4559",
				"QHEGZBAD EGZBADQH 02_2982_4267",
				"HQGBZDAE QGBZDAEH 03_3386_4575",
				"QBHDZGAE HDZGAEQB 02_2995_4312",
				"DQGZBEAH EAHDQGZB 04_3726_4878",
				"EGHDBZAQ DBZAQEGH 03_3360_4551",
				"DHEBZGAQ ZGAQDHEB 02_3021_4273",
				"BEGHZQAD HZQADBEG 02_2969_4242",
				"HBQEDGAZ GAZHBQED 03_3277_4474",
				"BQGHZDAE EBQGHZDA 03_3366_4553",
				"HBZQGEAD BZQGEADH 01_2611_4005",
				"EBGZQDAH BGZQDAHE 02_2928_4242",
				"EHBGZDAQ HBGZDAQE 02_2940_4241",
				"HZQDGBAE HZQDGBAE 02_3010_4278",
				"QZHBEGAD BEGADQZH 02_2989_4316",
				"BEHZGDAQ GDAQBEHZ 04_3711_4800",
				"EDGQHZAB GQHZABED 03_3378_4578",
				"HEDBQGAZ ZHEDBQGA 03_3396_4598",
				"HGZDQEAB HGZDQEAB 02_3002_4288",
				"HBZGEDAQ DAQHBZGE 03_3310_4487",
				"EQDHZGAB QDHZGABE 03_3341_4523",
				"HQZGDBAE AEHQZGDB 03_3291_4523",
				"GZBHEDAQ HEDAQGZB 02_2958_4218",
				"HEZDBQAG DBQAGHEZ 03_3335_4498",
				"DEGZQBAH ZQBAHDEG 05_4162_5170",
				"QZEHBGAD GADQZEHB 02_2978_4251",
				"QEZDHBAG AGQEZDHB 04_3720_4800",
				"QZDHBGAE BGAEQZDH 03_3320_4511",
				"GDZEQBAH AHGDZEQB 02_2987_4315",
				"DHZEGBAQ EGBAQDHZ 03_3327_4470",
				"HBDQEGAZ AZHBDQEG 03_3329_4514",
				"DHEZQBAG HEZQBAGD 03_3284_4477",
				"ZQDGEBAH GEBAHZQD 03_3441_4651",
				"EZDQHBAG EZDQHBAG 04_3710_4796",
				"ZDQGHBAE HBAEZDQG 05_4089_5107",
				"HZBGEDAQ QHZBGEDA 02_3038_4282",
				"ZHDBEGAQ HDBEGAQZ 03_3380_4571",
				"GZDBQEAH EAHGZDBQ 03_3373_4569",
				"DQHBGEAZ BGEAZDQH 04_3752_4823",
				"QBZHDEAG QBZHDEAG 04_3724_4782",
				"QBGDZEAH AHQBGDZE 02_3029_4337",
				"HZGBEDAQ EDAQHZGB 02_3014_4284",
				"DBQZHEAG AGDBQZHE 02_3015_4273",
				"GQHZDEAB GQHZDEAB 06_4498_5390",
				"DQEZHGAB EZHGABDQ 03_3372_4540",
				"EGHQDBAZ DBAZEGHQ 02_2967_4230",
				"EZGQDHAB ZGQDHABE 02_3072_4365",
				"HQBEZDAG AGHQBEZD 07_4791_5643",
				"HQEBZDAG EBZDAGHQ 04_3755_4825",
				"ZHDBGEAQ ZHDBGEAQ 04_3850_4930",
				"BQEZGDAH AHBQEZGD 03_3369_4563",
				"BEHGZDAQ AQBEHGZD 05_4054_5110",
				"BDEGQHAZ ZBDEGQHA 03_3345_4540",
				"BHQZDEAG DEAGBHQZ 03_3360_4536",
				"DBQGEHAZ QGEHAZDB 01_2641_3954",
				"BDGHQEAZ EAZBDGHQ 01_2620_4004",
				"GBQHDEAZ ZGBQHDEA 03_3440_4596",
				"HGBZEQAD HGBZEQAD 02_2972_4278",
				"HZEQDBAG BAGHZEQD 02_2963_4261",
				"BZHGQDAE HGQDAEBZ 03_3422_4603",
				"EZQDGBAH DGBAHEZQ 01_2698_4055",
				"QGBEZHAD GBEZHADQ 02_2965_4226",
				"EGZDHQAB GZDHQABE 01_2644_3997",
				"HDZQBEAG AGHDZQBE 04_3728_4862",
				"ZEGBHQAD EGBHQADZ 02_3031_4329",
				"GBEZQDAH EZQDAHGB 02_3002_4263",
				"ZDHQBEAG EAGZDHQB 05_4076_5126",
				"QHGZDBAE EQHGZDBA 04_3743_4846",
				"EGDBQHAZ BQHAZEGD 03_3377_4554",
				"ZBQHGDAE QHGDAEZB 03_3443_4630",
				"BZEGDHAQ DHAQBZEG 03_3349_4494",
				"GHQZEDAB DABGHQZE 02_2998_4273",
				"BDZGHEAQ DZGHEAQB 01_2647_4020",
				"QDBGEZAH BGEZAHQD 01_2643_4031",
				"HQGZBEAD DHQGZBEA 03_3386_4583",
				"EDBQHGAZ GAZEDBQH 02_2954_4280",
				"HDBQZGAE ZGAEHDBQ 04_3763_4886",
				"DGEZQBAH DGEZQBAH 02_2934_4216",
				"GQBZEHAD BZEHADGQ 03_3407_4554",
				"EQDHGBAZ ZEQDHGBA 03_3458_4640",
				"QHDEGZAB GZABQHDE 02_2980_4217",
				"EDGZBHAQ BHAQEDGZ 03_3342_4533",
				"ZGQEDHAB DHABZGQE 03_3373_4560",
				"EDZBQHAG EDZBQHAG 04_3662_4756",
				"BHQEGZAD BHQEGZAD 04_3664_4760",
				"EQBHGDAZ EQBHGDAZ 06_4476_5332",
				"QGHZEDAB BQGHZEDA 04_3754_4830",
				"EDQGHBAZ AZEDQGHB 01_2637_4001",
				"EQZGBHAD QZGBHADE 03_3372_4575",
				"QBZDHGAE GAEQBZDH 02_2981_4287",
				"BZQDGEAH HBZQDGEA 01_2641_4030",
				"ZEGQDBAH ZEGQDBAH 04_3764_4976",
				"DGQBHZAE QBHZAEDG 03_3414_4569",
				"BGDQHZAE ZAEBGDQH 03_3369_4577",
				"EQZHBGAD HBGADEQZ 03_3334_4545",
				"BZHQDGAE BZHQDGAE 04_3656_4804",
				"GEHQDBAZ AZGEHQDB 03_3321_4535",
				"DEQGBHAZ GBHAZDEQ 04_3682_4750",
				"GDQEBZAH ZAHGDQEB 05_4081_5112",
				"BHEGQZAD ZADBHEGQ 02_3051_4336",
				"HGEQZBAD ADHGEQZB 02_2920_4238",
				"ZQGDHEAB HEABZQGD 03_3484_4705",
				"DZBHEQAG EQAGDZBH 03_3440_4601",
				"DHZEBGAQ ZEBGAQDH 04_3747_4891",
				"EHDBGQAZ QAZEHDBG 04_3722_4794",
				"GDQHZEAB ZEABGDQH 03_3409_4596",
				"EZHBQGAD ADEZHBQG 03_3309_4513",
				"ZQHGBEAD ADZQHGBE 05_4145_5219",
				"HDGEZBAQ BAQHDGEZ 02_2941_4232",
				"GEHBQDAZ GEHBQDAZ 06_4410_5308",
				"QEDZGHAB EDZGHABQ 01_2599_3971",
				"GBZEDHAQ BZEDHAQG 03_3342_4507",
				"BHGDQZAE QZAEBHGD 03_3367_4570",
				"HGEQBZAD EQBZADHG 05_4097_5082",
				"BGZHQDAE GZHQDAEB 03_3297_4457",
				"EBZGDQAH HEBZGDQA 02_2924_4243",
				"ZGBQHDAE AEZGBQHD 04_3730_4879",
				"BEGHDZAQ EGHDZAQB 02_2959_4247",
				"DZBQGEAH DZBQGEAH 02_3060_4366",
				"DQBEGHAZ GHAZDQBE 03_3401_4563",
				"EBGQHDAZ ZEBGQHDA 03_3374_4597",
				"DZQHBGAE ZQHBGAED 05_4128_5119",
				"BQZDEHAG EHAGBQZD 04_3742_4847",
				"GQDBZHAE GQDBZHAE 04_3762_4846",
				"GBQZEDAH GBQZEDAH 04_3706_4798",
				"DGHBEZAQ ZAQDGHBE 03_3441_4614",
				"ZHQBEGAD DZHQBEGA 04_3800_4908",
				"DZEQBGAH HDZEQBGA 03_3296_4518",
				"HQDEBZAG QDEBZAGH 05_4090_5085",
				"BGQDZEAH QDZEAHBG 03_3431_4627",
				"ZEHBQGAD QGADZEHB 03_3433_4633",
				"GBZEHQAD GBZEHQAD 04_3652_4792",
				"GHBQDZAE HBQDZAEG 02_3011_4265",
				"BHDGEQAZ BHDGEQAZ 02_2996_4274",
				"ZDGHBEAQ GHBEAQZD 04_3785_4919",
				"ZEHGQBAD DZEHGQBA 03_3374_4568",
				"DHBQEGAZ QEGAZDHB 02_3012_4282",
				"QEZGHDAB ZGHDABQE 04_3770_4885",
				"DGHEQZAB GHEQZABD 03_3319_4494",
				"EGQDBHAZ ZEGQDBHA 03_3405_4613",
				"ZEGQHBAD BADZEGQH 02_2994_4362",
				"HBZDGEAQ GEAQHBZD 01_2630_4024",
				"QBGEZDAH ZDAHQBGE 04_3731_4858",
				"GBQZEHAD EHADGBQZ 01_2631_3980",
				"EGQHBZAD BZADEGQH 02_3004_4292",
				"ZDQBGEAH AHZDQBGE 02_3069_4376",
				"DBQHZEAG EAGDBQHZ 03_3387_4570",
				"HDEGZBAQ GZBAQHDE 01_2604_3957",
				"DHBEZGAQ DHBEZGAQ 04_3696_4804",
				"QBHZGEAD EADQBHZG 04_3654_4776",
				"DZGEHBAQ ZGEHBAQD 04_3727_4829",
				"EHBQDZAG AGEHBQDZ 03_3329_4490",
				"GEHQZDAB EHQZDABG 04_3693_4770",
				"BQDEZGAH GAHBQDEZ 05_4065_5077",
				"DZGEQHAB DZGEQHAB 04_3656_4772",
				"HZEGDQAB EGDQABHZ 01_2614_3967",
				"BQEGHZAD GHZADBQE 02_2998_4266",
				"DZBGHEAQ BGHEAQDZ 03_3409_4589",
				"ZDHEBGAQ QZDHEBGA 05_4105_5111",
				"HGDEQZAB ABHGDEQZ 03_3289_4492",
				"BHGDEZAQ AQBHGDEZ 03_3420_4574",
				"ZGBQEDAH EDAHZGBQ 03_3416_4619",
				"BGDEHQAZ DEHQAZBG 04_3739_4838",
				"QHEZBGAD QHEZBGAD 04_3656_4766",
				"DEGQBZAH AHDEGQBZ 03_3333_4568",
				"ZGBEDQAH EDQAHZGB 03_3405_4602",
				"HGDEZBAQ ZBAQHGDE 03_3346_4558",
				"QDHEZGAB ABQDHEZG 04_3726_4844",
				"QEDZBHAG DZBHAGQE 03_3350_4537",
				"GEQHDZAB HDZABGEQ 02_2941_4198",
				"DGQBEZAH BEZAHDGQ 02_3019_4278",
				"BHEGDQAZ HEGDQAZB 01_2599_3992",
				"BHZQDEAG GBHZQDEA 04_3694_4772",
				"DGBHZEAQ QDGBHZEA 01_2706_4052",
				"QDGEHZAB EHZABQDG 04_3698_4789",
				"EGHQZDAB EGHQZDAB 06_4348_5268",
				"GBEQZHAD HADGBEQZ 02_2936_4243",
				"HDQZBGAE ZBGAEHDQ 03_3355_4569",
				"ZDGEBHAQ AQZDGEBH 05_4173_5231",
				"GDBEQZAH GDBEQZAH 04_3702_4838",
				"GZEDHQAB QABGZEDH 01_2699_4071",
				"HGDZQBAE DZQBAEHG 02_3036_4282",
				"DZBEHGAQ AQDZBEHG 05_4131_5128",
				"DQHBGZAE AEDQHBGZ 03_3370_4561",
				"HGQEZBAD DHGQEZBA 02_2999_4291",
				"ZBGEQHAD GEQHADZB 04_3708_4832",
				"GDZBQEAH HGDZBQEA 02_2983_4276",
				"BEZDGQAH BEZDGQAH 02_2978_4328",
				"EZHDGBAQ HDGBAQEZ 02_3026_4301",
				"QBZEGHAD ADQBZEGH 02_2977_4288",
				"BHQEDGAZ AZBHQEDG 03_3362_4530",
				"QZEBGHAD DQZEBGHA 04_3722_4801",
				"DZGBHQAE EDZGBHQA 02_2982_4298",
				"HEBGQDAZ QDAZHEBG 03_3409_4587",
				"EHDZBGAQ AQEHDZBG 04_3715_4785",
				"BZEQGDAH ZEQGDAHB 03_3368_4557",
				"QEBZGHAD ZGHADQEB 03_3363_4546",
				"BQHDEZAG BQHDEZAG 08_5284_5982",
				"DZBGHQAE GHQAEDZB 02_3027_4278",
				"ZHBDGQAE EZHBDGQA 03_3417_4622",
				"QGBZHEAD HEADQGBZ 01_2698_4081",
				"EQDGBZAH EQDGBZAH 06_4450_5402",
				"EDQHZBAG QHZBAGED 02_3000_4257",
				"DBHEGQAZ ZDBHEGQA 04_3730_4841",
				"DQZBGEAH ZBGEAHDQ 04_3803_4921",
				"BHEZGDAQ QBHEZGDA 04_3654_4735",
				"HZDEGQAB EGQABHZD 03_3277_4467",
				"EQHGDZAB ZABEQHGD 05_4093_5122",
				"DEHGQBAZ BAZDEHGQ 02_2976_4265",
				"BZGHEQAD GHEQADBZ 02_3037_4303",
				"BZEHQDAG BZEHQDAG 06_4398_5242",
				"ZEQDGHAB HABZEQDG 03_3394_4582",
				"BZDGHQAE AEBZDGHQ 01_2618_4019",
				"GEDQHBAZ QHBAZGED 01_2639_4004",
				"DQBZGEAH QBZGEAHD 02_3042_4323",
				"EZGHBQAD EZGHBQAD 04_3666_4810",
				"EZDGHBAQ BAQEZDGH 03_3331_4533",
				"HDEGBQAZ QAZHDEGB 02_2967_4242",
				"EQDGZHAB ABEQDGZH 02_2933_4256",
				"BGDQEZAH HBGDQEZA 01_2645_4062",
				"ZHDGBQAE DGBQAEZH 03_3429_4666",
				"ZGEHDQAB BZGEHDQA 04_3734_4831",
				"EZDGHQAB HQABEZDG 03_3403_4554",
				"QEGZDHAB ABQEGZDH 02_2932_4238",
				"GDQHBZAE AEGDQHBZ 02_3006_4304",
				"EDGHBQAZ DGHBQAZE 03_3312_4518",
				"BZGDEHAQ HAQBZGDE 01_2644_4032",
				"QDZGEBAH QDZGEBAH 02_2990_4310",
				"ZHEGBQAD BQADZHEG 05_4160_5156",
				"DGQHBZAE DGQHBZAE 04_3690_4754",
				"QDGBEZAH DGBEZAHQ 02_3032_4311",
				"QBEHZDAG BEHZDAGQ 04_3668_4732",
				"GDEQHBAZ DEQHBAZG 03_3336_4503",
				"HQEZDGAB DGABHQEZ 02_3018_4258",
				"DZBGEQAH QAHDZBGE 02_3038_4339",
				"DQBEHZAG ZAGDQBEH 05_4186_5208",
				"ZHGEQDAB BZHGEQDA 05_4097_5117",
				"HQBZEGAD GADHQBZE 02_2966_4265",
				"DEZGBHAQ DEZGBHAQ 04_3636_4764",
				"DHBGQEAZ AZDHBGQE 02_2986_4277",
				"QGDHZBAE HZBAEQGD 01_2635_3989"};		
		
		for (int c=0 ; c<runs.length ; c++) {
			String[] passTo = new String[3];
			passTo[0] = runs[c].substring(0, 8);
			passTo[1] = runs[c].substring(9, 17);
			passTo[2] = runs[c].substring(18);
			maintmp(passTo);
		}

		
		
		
		
	}
	
	public static void maintmp(String[] args) throws Exception {
		init(args);
		
		// Defining the arrangement 
		String topTrue = args[0]; 
//		String topTrue = "ZQHEBDAG"; // OMS
//		String topTrue = "ZQHEBDAG"; // OMS+1
//		String topTrue = "ZQHEBDAG"; // OMS-1
//		String topTrue = "ZEAHDQGB"; // Willison
//		String topTrue = "AZBGQDEH"; // Yao
		String chainsTop = "FEAGDHCB";
		String chainsTop_E = "feagdhcb";
		String genesTop =  "ZEAHDQGB";

		String botTrue = args[1]; 
		String botToWrite = ""+botTrue.charAt(0);
		for (int c=7 ; c>0 ; c--) {
			botToWrite = botToWrite + botTrue.charAt(c);
		}
//		String botTrue = "ZQHEBDAG"; // OMS
//		String botTrue = "GZQHEBDA"; // OMS+1
//		String botTrue = "QHEBDAGZ"; // OMS-1
//		String botTrue = "ZEAHDQGB"; // Willison
//		String botTrue = "AZBGQDEH"; // Yao
		String chainsBot = "NMIOLPKJ";
		String chainsBot_E = "nmiolpkj";
		String genesBot =  "ZEAHDQGB";

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
		int modelNumber=36; //int modelNumber = 0;
		for (int delta1=4 ; delta1<5 ; delta1++) {//(int delta1=0 ; delta1<8 ; delta1++) { // moving the real 2-fold "delta1" units in the top ring order from willisons 2-fold. In complex 3P9E
			for (int delta2=4 ; delta2<5 ; delta2++) {//(int delta2=0 ; delta2<8 ; delta2++) { // moving the real 2-fold "delta1" units in the top ring order from willisons 2-fold. In complex 3P9E
				String outputFile;
				outputFile = "chain2gene_"+topTrue+"_"+botToWrite+"_"+args[2]+".txt";
//				if (modelNumber>9) {
//					outputFile = "chain2gene_"+modelNumber+".txt";
//				}
//				else {
//					outputFile = "chain2gene_0"+modelNumber+".txt";
//				}
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
		modelNumber=36; //int modelNumber = 0;
		for (int delta1=4 ; delta1<5 ; delta1++) {//(int delta1=0 ; delta1<8 ; delta1++) { // moving the real 2-fold "delta1" units in the top ring order from willisons 2-fold. In complex 3P9E
			for (int delta2=4 ; delta2<5 ; delta2++) {//(int delta2=0 ; delta2<8 ; delta2++) { // moving the real 2-fold "delta1" units in the top ring order from willisons 2-fold. In complex 3P9E
					String outputFile;
					outputFile = "model_"+topTrue+"_"+botToWrite+"_"+args[2]+".pdb";
//					if (modelNumber>9) {
//						outputFile = "model_"+modelNumber+"_"+topTrue+"_"+botTrue+"_"+delta1+"_"+delta2+".pdb";
//					}
//					else {
//						outputFile = "model_0"+modelNumber+"_"+topTrue+"_"+botTrue+"_"+delta1+"_"+delta2+".pdb";
//					}
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
		AtomList newListSCWRL = new AtomList("/Users/nirka/TRiC/New_Work_on_TRiC_14_8_2011/Gunnar/HM_Yeast_"+queryUnit+".scwrl.pdb");
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
