#procedure tadpole(TOPO,i)
if (count(`TOPO',1));
#switch `TOPO'
#case Md6
#case Md15d
#case Md17d
if (count(`i'1,1) >= 0) Discard;
#break
#case Md11
#case Md15M2
#case Md34d
#case Md14M6
#case Md28d
if (((count(`i'2,1) >= 0))||((count(`i'1,1) >= 0))) Discard;
#break
#case Md15
#case Mt1star05M1
#case Mt1star05M2
#case Mt1star05M3
#case Mt1star05M4
#case Md110d
#case Md30M7
#case Md33d
#case Md77d
#case Md78d
#case Md85d
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))) Discard;
#break
#case Md15c
if (count(`i'1,1) >= 0) Discard;
#break
#case Md14
#case Md30M2
#case Md50d
#case Md87d
#case Md111d
#case Md107d
if (((count(`i'2,1) >= 0))||((count(`i'1,1) >= 0))) Discard;
#break
#case Md28
#case Md33M2
#case Md50M2
if (((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md28c
if (((count(`i'3,1) >= 0))||((count(`i'2,1) >= 0))) Discard;
#break
#case Md50
#case Md81M1
#case Md77M2
#case Md87M2
#case Md87M3
#case Md60M2
#case Md60M4
#case Md60M5
#case Md60M8
#case Md78M5
#case Md107M1
#case Md107M5
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))) Discard;
#break
#case Md50c
if (((count(`i'4,1) >= 0))||((count(`i'3,1) >= 0))) Discard;
#break
#case Md33
#case Md92M1
#case Md92M2
#case Md92M3
#case Md92M4
#case Md81M3
#case Md77M3
#case Md172d
#case Md78M4
#case Md188d
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))) Discard;
#break
#case Md33c
if (((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md77
#case Md137M1
#case Md137M2
#case Md118M1
#case Md118M6
#case Md221d
#case Md121M1
#case Md121M5
#case Md190M2
#case Md190M5
#case Md247d
#case Md148M1
if (((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md77c
if (((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md4
#case Md25d
#case Md9d
#case Md22d
#case Md13d
if (count(`i'1,1) >= 0) Discard;
#break
#case Md9
#case Md42d
#case Md22M2
if (((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'1,1) >= 0))) Discard;
#break
#case Md9c
if (count(`i'2,1) >= 0) Discard;
#break
#case Md22
#case Md75d
#case Md43M1
#case Md43M3
#case Md43M4
#case Md43M7
#case Md79d
#case Md55d
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md22c
if (count(`i'2,1) >= 0) Discard;
#break
#case Md3
#case Md24d
#case Md23d
#case Md26d
if (((count(`i'2,1) >= 0))||((count(`i'1,1) >= 0))) Discard;
#break
#case Md26
#case Md43M2
#case Md54d
if (((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md26c
if (((count(`i'2,1) >= 0))||((count(`i'1,1) >= 0))) Discard;
#break
#case Md43
#case Md72d
#case Md60d
#case Md124d
#case Md121d
#case Md150d
#case Md126d
#case Md149d
if (((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md17
#case Md30M1
#case Md37d
#case Md91d
#case Md80d
#case Md106d
if (((count(`i'2,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md17c
if (count(`i'2,1) >= 0) Discard;
#break
#case Md37
#case Md81M4
#case Md60M3
#case Md180d
#case Md80M5
#case Md106M5
if (((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md37c
if (((count(`i'4,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))) Discard;
#break
#case Md60
#case Md118M2
#case Md118M3
#case Md281d
#case Md121M6
#case Md244d
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md60c
if (((count(`i'4,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'3,1) >= 0))) Discard;
#break
#case Md13
#case Md30M3
#case Md32d
#case Md112d
#case Md105d
if (((count(`i'2,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md13c
if (count(`i'2,1) >= 0) Discard;
#break
#case Md30
#case Md57d
#case Md81d
#case Md119d
#case Md140d
#case Md191d
#case Md134d
#case Md117d
#case Md148d
if (((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md32
#case Md79M5
#case Md81M8
#case Md189d
#case Md151d
#case Md55M5
if (((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md32c
if (((count(`i'3,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md81
#case Md137M3
#case Md137M5
#case Md118M4
#case Md118M9
#case Md235d
#case Md191M5
#case Md241d
#case Md148M6
if (((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md81c
if (((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'3,1) >= 0))) Discard;
#break
#case Md87
#case Md118M5
#case Md218d
#case Md191M2
#case Md233d
if (((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md87c
if (((count(`i'4,1) >= 0))||((count(`i'3,1) >= 0))) Discard;
#break
#case Md118
#case Mnostar0M1
#case Mnostar0M3
#case Mnostar0M4
#case Mnostar0M6
#case Md209d
#case Md306d
#case Md197M2
#case Md197M6
#case Md337d
#case Md366d
#case Md370d
if (((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md5
#case Md20d
#case Md16d
#case Md19d
#case Md21d
if (count(`i'1,1) >= 0) Discard;
#break
#case Md8
#case Md20M2
#case Md52d
#case Md27d
#case Md38d
#case Md29M1
if (((count(`i'2,1) >= 0))||((count(`i'1,1) >= 0))) Discard;
#break
#case Md20
#case Mt1star55M1
#case Mt1star55M2
#case Mt1star55M3
#case Mt1star55M4
#case Md83d
#case Md51d
#case Md89d
#case Md47M1
#case Md63d
#case Md61d
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))) Discard;
#break
#case Md20c
if (count(`i'1,1) >= 0) Discard;
#break
#case Md27
#case Md52M2
#case Md51M2
#case Md35M3
#case Md41M7
#case Md53d
if (((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md27c
if (((count(`i'3,1) >= 0))||((count(`i'2,1) >= 0))) Discard;
#break
#case Md51
#case Md83M2
#case Md93M4
#case Md104M3
#case Md86M1
#case Md66M1
#case Md66M3
#case Md66M4
#case Md66M8
#case Md94M8
#case Md102d
#case Md59M5
#case Md61M4
if (((count(`i'5,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md51c
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md52
#case Md83M3
#case Md76M4
#case Md57M8
#case Md58M2
#case Md58M3
#case Md58M4
#case Md58M8
#case Md104M5
#case Md103M5
#case Md61M3
#case Md84M5
if (((count(`i'5,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md52c
if (((count(`i'3,1) >= 0))||((count(`i'2,1) >= 0))) Discard;
#break
#case Md83
#case Mfastar3M5
#case Mbustar5M3
#case Md136M1
#case Md136M6
#case Md128M1
#case Md155d
#case Md140M3
#case Md215d
#case Md239d
#case Md145M2
#case Md145M4
#case Md147M1
#case Md147M4
#case Md116M6
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md83c
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md10
#case Md18M1
#case Md19M2
#case Md41d
if (count(`i'1,1) >= 0) Discard;
#break
#case Md18
#case Mt1star45M3
#case Mt1star24M1
#case Mt1star24M3
#case Md35d
#case Md103d
#case Md73d
#case Md62d
#case Md84d
if (count(`i'2,1) >= 0) Discard;
#break
#case Md41
#case Md57M2
#case Md103M8
#case Md73M1
#case Md86M3
#case Md94M5
#case Md72M1
#case Md72M3
#case Md72M5
#case Md72M8
#case Md100d
#case Md159d
#case Md173d
#case Md59M2
#case Md62M8
#case Md84M2
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))) Discard;
#break
#case Md41c
if (count(`i'2,1) >= 0) Discard;
#break
#case Md35
#case Md76M2
#case Md93M2
#case Md73M2
#case Md96d
#case Md62M5
#case Md183d
if (((count(`i'5,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'2,1) >= 0))) Discard;
#break
#case Md35c
if (count(`i'4,1) >= 0) Discard;
#break
#case Md73
#case Md144M4
#case Md127M9
#case Md163d
#case Md216d
#case Md146M1
#case Md124M2
#case Md124M6
#case Md122M2
#case Md255d
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'3,1) >= 0))) Discard;
#break
#case Md73c
if (count(`i'4,1) >= 0) Discard;
#break
#case Md53
#case Md96M3
#case Md102M3
#case Md100M8
if (((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md53c
if (((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md100
#case Md163M3
#case Md170M3
#case Md160M2
#case Md159M2
#case Md159M3
#case Md166M2
#case Md166M4
#case Md166M5
#case Md166M6
#case Md167M3
#case Md171M2
#case Md173M2
#case Md173M4
#case Md165M2
#case Md162M6
#case Md183M2
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))) Discard;
#break
#case Md100c
if (((count(`i'3,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md96
#case Md176M4
#case Md163M9
#case Md160M5
#case Md171M5
#case Md155M2
#case Md187M1
#case Md187M3
#case Md162M5
#case Md164M1
#case Md164M2
#case Md164M6
#case Md164M9
#case Md154M3
#case Md183M5
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md96c
if (((count(`i'5,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'3,1) >= 0))) Discard;
#break
#case Md163
#case Md213M4
#case Md209M4
#case Md209M5
#case Md215M2
#case Md216M1
#case Md263M1
#case Md284M2
#case Md284M6
#case Md257M2
#case Md257M6
#case Md211M1
#case Md211M2
#case Md219M2
#case Md266M10
#case Md279M2
#case Md248M2
#case Md255M2
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md163c
if (((count(`i'3,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md102
#case Md176M3
#case Md170M9
#case Md167M5
#case Md156M1
#case Md156M2
#case Md156M3
#case Md156M4
#case Md155M3
#case Md165M6
#case Md154M4
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md102c
if (((count(`i'5,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'3,1) >= 0))) Discard;
#break
#case Md155
#case Md207M1
#case Md217M1
#case Md217M7
#case Md220M1
#case Md280M4
#case Md215M5
#case Md239M6
#case Md205M10
#case Md287M2
#case Md287M4
#case Md282M1
#case Md282M4
#case Md225M6
if (((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md155c
if (((count(`i'6,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'1,1) >= 0))) Discard;
#break
#case Md19
#case Mt1star45M1
#case Md94d
#case Md59d
if (count(`i'1,1) >= 0) Discard;
#break
#case Md19c
if (count(`i'1,1) >= 0) Discard;
#break
#case Md94
#case Md128M4
#case Md135M6
#case Md167d
#case Md212d
#case Md141M2
#case Md134M4
#case Md123M2
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))) Discard;
#break
#case Md94c
if (count(`i'1,1) >= 0) Discard;
#break
#case Md167
#case Md206M7
#case Md212M7
#case Md212M10
#case MbebeM2
#case MbebeM4
#case Md207M2
#case Md234M2
#case Md270M3
#case Md270M7
#case Md281M4
#case Md281M6
#case Md215M10
#case Md264M7
#case Md254M2
#case Md278M4
#case Md223M2
#case Md277M6
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md167c
if (((count(`i'5,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md215
#case Md378M1
#case Md378M4
#case Md306M1
#case Md306M5
#case Md381M4
#case Md381M11
#case Md321M3
#case Md321M7
#case Md305M3
#case Md345M11
#case Md310M1
#case Md334M2
#case Md333M2
#case Md338M4
#case Md323M2
#case Md323M4
#case Md299M3
#case Md300M2
#case Md342M3
#case Md342M4
#case Md317M7
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md215c
if (((count(`i'4,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))) Discard;
#break
#case Md12
#case Md49d
#case Md45M2
#case Md45M3
#case Md64d
#case Md113d
#case Md90d
#case Md108d
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md49
#case Md98M2
#case Md72M2
#case Md113M2
#case Md99d
#case Md108M8
if (((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md49c
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md72
#case Md135M5
#case Md166d
#case Md124M9
#case Md256d
if (((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md72c
if (((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md99
#case Md160M3
#case Md166M9
#case Md185M3
#case Md185M6
#case Md157M2
#case Md158M5
#case Md174M9
#case Md179M5
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md99c
if (((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md166
#case Md209M2
#case Md209M3
#case Md206M3
#case Md281M2
#case Md257M7
#case Md256M6
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md166c
if (((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'2,1) >= 0))) Discard;
#break
#case Md38
#case Md89M1
#case Md98M3
#case Md57M1
#case Md64M1
#case Md90M1
#case Md63M2
#case Md158d
#case Md179d
if (((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md38c
if (((count(`i'3,1) >= 0))||((count(`i'2,1) >= 0))) Discard;
#break
#case Md63
#case Md157d
#case Md140M4
#case Md262d
#case Md131M9
#case Md274d
#case Md138M1
#case Md138M4
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md63c
if (((count(`i'3,1) >= 0))||((count(`i'1,1) >= 0))) Discard;
#break
#case Md157
#case Md208M1
#case Md208M4
#case Md280M5
#case Md281M5
#case Md205M4
#case Md229M7
#case Md262M6
#case Md278M3
#case Md226M10
#case Md274M6
#case Md286M1
#case Md286M2
#case Md286M5
#case Md286M10
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md157c
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md59
#case Md165d
#case Md175M1
#case Md140M6
#case Md146M9
#case Md264d
#case Md283d
if (((count(`i'6,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'3,1) >= 0))) Discard;
#break
#case Md59c
if (count(`i'3,1) >= 0) Discard;
#break
#case Md165
#case Md271M1
#case Md270M1
#case Md270M4
#case Md280M7
#case Md281M7
#case Md281M10
#case Md239M7
#case Md236M10
#case Md263M10
#case Md264M2
#case Md279M3
#case Md283M3
#case Md283M4
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md165c
if (((count(`i'6,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'3,1) >= 0))) Discard;
#break
#case Md281
#case Md308M1
#case Md306M2
#case Md306M11
#case Md321M8
#case Md305M11
#case Md337M3
#case Md337M4
#case Md341M8
#case Md322M7
#case Md307M1
#case Md307M2
#case Md307M5
#case Md307M11
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md281c
if (((count(`i'1,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md159
#case Md209M6
#case Md212M6
#case McrossM5
#case McrossM6
#case Md216M2
#case Md236M1
#case Md285M2
#case Md256M2
#case Md256M10
#case Md283M2
#case Md238M4
if (((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md159c
if (((count(`i'4,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'5,1) >= 0))) Discard;
#break
#case Md16
#case Mt1star45M2
#case Mt1star34M1
#case Mt1star34M2
#case Md69d
#case Md36d
#case Md39d
#case Md68d
#case Md109d
#case Md82d
if (((count(`i'2,1) >= 0))||((count(`i'1,1) >= 0))) Discard;
#break
#case Md16c
if (count(`i'1,1) >= 0) Discard;
#break
#case Md39
#case Md57M4
#case Md104M8
#case Md71M1
#case Md86M2
#case Md178d
#case Md65M5
#case Md109M2
#case Md82M2
if (((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'1,1) >= 0))) Discard;
#break
#case Md39c
if (((count(`i'4,1) >= 0))||((count(`i'2,1) >= 0))) Discard;
#break
#case Md57
#case Mbustar5M9
#case Md127M1
#case Md136M2
#case Md136M3
#case Md135M1
#case Md160d
#case Md140M5
#case Md285d
#case Md266d
#case Md134M2
if (((count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md57c
if (((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md36
#case Md69M2
#case Md93M3
#case Md71M2
#case Md101d
#case Md68M5
#case Md169d
#case Md182d
if (((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md36c
if (((count(`i'4,1) >= 0))||((count(`i'3,1) >= 0))) Discard;
#break
#case Md101
#case Md176M6
#case Md181M3
#case Md181M9
#case Md170M5
#case Md160M9
#case Md178M3
#case Md168M6
#case Md161M5
#case Md169M2
#case Md182M2
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md101c
if (((count(`i'3,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md160
#case Md209M7
#case Md209M10
#case Md206M6
#case MbebeM1
#case MbebeM3
#case McrossM1
#case McrossM2
#case McrossM3
#case McrossM4
#case Md217M2
#case Md217M3
#case Md280M6
#case Md285M6
#case Md205M2
#case Md211M6
#case Md211M7
#case Md219M10
#case Md266M6
#case Md278M2
if (((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md160c
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md209
#case MbubuM4
#case Md306M3
#case Md381M1
#case Md381M2
#case Md305M1
#case Md336M2
#case Md336M7
#case Md339M2
#case Md339M11
#case Md337M2
#case Md384M1
#case Md384M2
#case Md384M4
#case Md384M5
#case Md384M6
#case Md384M7
#case Md383M2
#case Md383M3
#case Md383M4
#case Md383M5
if (((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md209c
if (((count(`i'7,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md21
#case Md47M2
#case Md88d
#case Md65d
#case Md67d
if (((count(`i'2,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md21c
if (count(`i'1,1) >= 0) Discard;
#break
#case Md65
#case Md168d
#case Md175M4
#case Md140M9
#case Md260d
#case Md130M6
#case Md272d
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md65c
if (((count(`i'2,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md140
#case Md280d
#case Md333d
#case Md199M10
#case Md318d
#case Md195M1
#case Md195M4
if (((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md140c
if (((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'2,1) >= 0))) Discard;
#break
#case Md168
#case Md271M4
#case Md280M10
#case Md260M6
#case Md259M7
#case Md272M1
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md168c
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md280
#case Md308M4
#case Md306M6
#case Md306M8
#case Md333M7
#case Md324M8
#case Md338M1
#case Md343M11
#case Md318M7
#case Md309M1
#case Md309M5
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md280c
if (((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md212
#case MbubuM1
#case Md306M7
#case Md344M3
#case Md310M6
#case Md363M4
#case Md312M7
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))) Discard;
#break
#case Md212c
if (((count(`i'2,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))) Discard;
#break
#case Md306
#case Mno6M1
#case Mno6M5
#case Md394M1
#case Md393M1
#case Md393M4
#case Md417M8
#case Md416M12
#case Md424M2
#case Md424M12
if (((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md306c
if (((count(`i'6,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Mt1star05
#case Md58d
#case Md92d
#case Md190d
#case Md147d
#case Md120d
if (((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md92
#case Md137M4
#case Md190M6
#case Md224d
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md92c
if (((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md79
#case Md137M6
#case Md137M9
#case Md191M6
#case Md242d
#case Md245d
#case Md149M2
#case Md149M6
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md79c
if (((count(`i'2,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md137
#case Mnostar0M2
#case Mnostar0M5
#case Mnostar0M7
#case Mnostar0M10
#case Md211d
#case Md197M7
#case Md338d
#case Md296d
#case Md192M2
#case Md192M7
if (((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Mnostar0
#case Md383d
#case Md393d
#case Md400d
if (((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md68
#case Md161d
#case Md146M2
#case Md130M9
#case Md126M2
#case Md126M4
#case Md268d
if (((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md68c
if (((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md161
#case Md236M2
#case Md285M7
#case Md263M2
#case Md211M3
#case Md211M4
#case Md259M10
#case Md272M3
#case Md267M3
#case Md240M4
#case Md265M2
#case Md265M4
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md161c
if (((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md58
#case Mbustar5M6
#case Md164d
#case Md284d
#case Md147M3
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md58c
if (((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))) Discard;
#break
#case Md164
#case Md205M1
#case Md284M7
#case Md211M5
#case Md282M3
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md164c
if (((count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'2,1) >= 0))) Discard;
#break
#case Md211
#case Md305M2
#case Md336M8
#case Md324M11
#case Md338M3
#case Md383M6
#case Md383M7
#case Md383M8
#case Md383M11
#case Md340M2
#case Md340M7
if (((count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md211c
if (((count(`i'5,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md383
#case Md394M2
#case Md393M3
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md383c
if (((count(`i'3,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))) Discard;
#break
#case Md25
#case Mt1star34M7
#case Mt1star24M7
#case Md56d
if (((count(`i'2,1) >= 0))||((count(`i'1,1) >= 0))) Discard;
#break
#case Md25c
if (count(`i'1,1) >= 0) Discard;
#break
#case Mt1star34
#case Md71d
#case Md130d
#case Md125d
if (((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md42
#case Md75M3
#case Md76M5
#case Md71M5
#case Md97d
#case Md152d
#case Md56M5
if (((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md42c
if (((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'1,1) >= 0))) Discard;
#break
#case Md71
#case Md144M3
#case Md127M4
#case Md181d
#case Md130M5
#case Md258d
#case Md125M2
if (((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md71c
if (((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md97
#case Md181M5
#case Md153M5
#case Md187M6
#case Md152M6
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md97c
if (((count(`i'3,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'2,1) >= 0))) Discard;
#break
#case Md181
#case Md213M6
#case McrossM7
#case McrossM10
#case Md219M6
#case Md259M5
#case Md258M2
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md181c
if (((count(`i'3,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Mcross
#case MbubuM7
#case Md339M7
#case Md324M6
#case Md384M3
#case Md384M8
#case Md384M11
if (((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md384
#case Md394M3
if (((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md23
#case Md45M7
#case Md74d
if (((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))) Discard;
#break
#case Md23c
if (((count(`i'2,1) >= 0))||((count(`i'1,1) >= 0))) Discard;
#break
#case Md45
#case Md98d
#case Md115d
#case Md143d
if (((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md34
#case Md110M3
#case Md98M8
#case Md74M1
#case Md86M8
#case Md111M2
#case Md95d
#case Md186d
#case Md85M3
if (((count(`i'5,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md34c
if (((count(`i'3,1) >= 0))||((count(`i'2,1) >= 0))) Discard;
#break
#case Md98
#case Mfastar2M1
#case Mfastar2M3
#case Md135M2
#case Md115M2
#case Md185d
#case Md143M6
if (((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md98c
if (((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md103
#case Mbustar5M2
#case Md135M3
#case Md171d
#case Md175M3
#case Md142M6
#case Md133M4
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))) Discard;
#break
#case Md103c
if (count(`i'1,1) >= 0) Discard;
#break
#case Md29
#case Md47M3
if (((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md24
#case Md70d
#case Md47M7
if (((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md24c
if (((count(`i'2,1) >= 0))||((count(`i'1,1) >= 0))) Discard;
#break
#case Md47
#case Md86d
#case Md114d
#case Md133d
#case Md116d
#case Md131d
if (((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md86
#case Md144M1
#case Md136M4
#case Md136M9
#case Md128M6
#case Md135M9
#case Md170d
#case Md234d
#case Md236d
#case Md133M2
#case Md116M3
if (((count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md86c
if (((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md135
#case Mnostar5M4
#case Mnostar5M6
#case Md206d
#case Md344d
#case Md332d
#case Md196M7
#case Md198M2
#case Md198M4
if (((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md171
#case Md206M2
#case Md271M3
#case Md216M7
#case Md205M7
#case Md250M6
#case Md275M6
#case Md276M4
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md171c
if (((count(`i'5,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md95
#case Md184M3
#case Md170M6
#case Md185M5
#case Md172M5
#case Md186M2
#case Md188M3
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md95c
if (((count(`i'5,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md170
#case Md213M3
#case Md206M5
#case Md207M5
#case Md217M4
#case Md217M10
#case Md234M4
#case Md236M6
#case Md276M2
#case Md225M3
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md170c
if (((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md185
#case Md206M10
#case MbebeM5
#case Md208M3
#case Md208M6
#case Md232M6
#case Md252M6
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md185c
if (((count(`i'2,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md206
#case MbubuM2
#case MbubuM3
#case Md308M3
#case Md381M3
#case Md381M8
#case Md305M7
#case Md344M7
#case Md385M2
#case Md385M7
#case Md332M7
#case Md328M7
#case Md347M2
#case Md347M4
if (((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md206c
if (((count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Mt1star45
#case Md93d
#case Md104d
#case Md175d
#case Md146d
#case Md141d
#case Md123d
if (((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))) Discard;
#break
#case Md93
#case Mfastar3M2
#case Md144M2
#case Md176d
#case Md146M6
#case Md141M6
#case Md223d
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md93c
if (((count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'2,1) >= 0))) Discard;
#break
#case Md75
#case Mfastar3M6
#case Md144M6
#case Md153d
#case Md222d
#case Md142M2
#case Md130M2
#case Md150M3
#case Md150M6
#case Md249d
if (((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md75c
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md74
#case Md144M9
#case Md184d
#case Md143M2
#case Md117M3
#case Md251d
if (((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md74c
if (((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))) Discard;
#break
#case Md144
#case Mnostar6M5
#case Mnostar6M10
#case Md213d
#case Md325d
#case Md196M2
#case Md193M3
#case Md303d
if (((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))) Discard;
#break
#case Md153
#case Md213M5
#case Md222M1
#case Md222M4
#case Md220M7
#case Md250M2
#case Md259M2
#case Md289M3
#case Md289M7
#case Md249M6
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md153c
if (((count(`i'5,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md184
#case Md213M7
#case MbebeM10
#case Md218M3
#case Md218M4
#case Md221M6
#case Md236M7
#case Md244M3
#case Md244M7
#case Md252M2
#case Md241M3
#case Md251M2
#case Md237M4
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md184c
if (((count(`i'4,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md176
#case Md213M10
#case MbebeM6
#case MbebeM7
#case Md207M3
#case Md220M4
#case Md271M7
#case Md205M3
#case Md263M7
#case Md254M7
#case Md223M7
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md176c
if (((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'2,1) >= 0))) Discard;
#break
#case Md213
#case MbubuM5
#case MbubuM6
#case Md378M3
#case Md378M7
#case Md379M3
#case Md379M7
#case Md345M3
#case Md325M6
#case Md328M2
#case Md324M2
#case Md302M3
#case Md303M2
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md213c
if (((count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Mbebe
#case MbubuM8
#case MbubuM11
#case MnonoM3
#case MnonoM4
#case Md308M7
#case Md305M8
#case Md345M7
if (((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Mbubu
#case Mno6M4
#case Mno6M8
#case Md391M3
#case Md391M8
#case Md387M7
#case Md394M4
if (((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md216
#case Md381M5
#case Md381M7
#case Md311M7
#case Md344M11
#case Md325M1
#case Md327M2
#case Md324M1
#case Md326M1
#case Md326M7
#case Md314M2
#case Md320M2
#case Md365M5
if (((count(`i'8,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md216c
if (((count(`i'6,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'3,1) >= 0))) Discard;
#break
#case Md136
#case Mnostar6M1
#case Md217d
#case Md316d
#case Md194M2
#case Md342d
if (((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md217
#case MnonoM1
#case MnonoM6
#case Md381M6
#case Md316M7
#case Md379M11
#case Md324M7
#case Md315M2
#case Md342M7
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md217c
if (((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md381
#case Md392M8
#case Md391M12
#case Md387M1
#case Md386M2
#case Md396M2
#case Md394M5
#case Md394M6
if (((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md146
#case Md263d
#case Md320d
#case Md199M3
if (((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md146c
if (((count(`i'5,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))) Discard;
#break
#case Md263
#case Md311M3
#case Md305M4
#case Md310M11
#case Md320M4
#case Md343M3
#case Md359M4
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))) Discard;
#break
#case Md263c
if (((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md236
#case Md316M3
#case Md316M8
#case Md305M5
#case Md325M11
#case Md310M2
#case Md332M2
#case Md324M3
#case Md314M4
#case Md372M4
#case Md341M2
#case Md341M4
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md236c
if (((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md89
#case Mfastar2M4
#case Mfastar2M5
#case Mbustar5M1
#case Md114M1
#case Md229d
#case Md134M3
#case Md138M2
#case Md138M9
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md89c
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md104
#case Mbustar5M4
#case Md128M3
#case Md175M6
#case Md123M6
if (((count(`i'2,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md104c
if (((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Mbustar5
#case Mnostar5M1
#case Mnostar5M3
#case Md311d
#case Md205d
#case Md334d
#case Md194M4
#case Md195M2
#case Md195M10
if (((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md205
#case Md311M8
#case Md305M6
#case Md385M1
#case Md385M11
#case Md334M7
#case Md324M5
#case Md315M4
#case Md309M2
#case Md309M11
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))) Discard;
#break
#case Md205c
if (((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))) Discard;
#break
#case Md305
#case Md392M3
#case Md392M9
#case Md387M12
#case Md396M4
#case Md394M7
#case Md394M8
#case Md416M3
#case Md424M1
#case Md424M6
if (((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md130
#case Md259d
#case Md352d
#case Md199M5
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md130c
if (((count(`i'5,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'2,1) >= 0))) Discard;
#break
#case Md259
#case Md325M2
#case Md327M7
#case Md324M4
#case Md352M2
#case Md352M4
#case Md361M5
#case Md343M5
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md259c
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md324
#case Md387M2
#case Md386M8
#case Md397M2
#case Md397M4
#case Md394M9
#case Md394M12
#case Md416M8
#case Md420M1
#case Md420M2
#case Md420M6
#case Md420M12
if (((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md394
#case Mno3M2
#case Mno3M4
#case Mno2M1
#case Mno2M4
#case Mno2M5
#case Mno2M8
if (((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))) Discard;
#break
#case Md172
#case Md207M7
#case Md235M2
#case Md208M2
#case Md208M5
#case Md221M5
#case Md246M7
#case Md276M5
#case Md224M3
#case Md224M7
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md172c
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))) Discard;
#break
#case Md221
#case Md308M2
#case Md378M2
#case Md378M8
#case Md349M3
#case Md310M3
#case Md356M3
#case Md356M8
#case Md314M5
#case Md295M3
#case Md296M3
#case Md370M3
#case Md370M4
#case Md368M4
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md221c
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md222
#case MlalaM8
#case Md378M5
#case Md353M3
#case Md325M3
#case Md352M5
#case Md293M3
#case Md293M7
#case Md330M6
if (((count(`i'2,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md222c
if (((count(`i'4,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md110
#case Mfastar2M6
#case Mfastar2M9
#case Md128M2
#case Md119M2
#case Md133M5
#case Md120M3
#case Md120M6
if (((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md110c
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Mt1star55
#case Md66d
#case Md145d
#case Md138d
if (((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md66
#case Md128M9
#case Md156d
#case Md270d
#case Md145M5
if (((count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md66c
if (((count(`i'5,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'1,1) >= 0))) Discard;
#break
#case Md128
#case Mnostar5M7
#case Mnostar5M10
#case Md207d
#case Md354d
#case Md310d
#case Md194M5
#case Md201M3
#case Md201M7
if (((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md156
#case Md207M4
#case Md270M5
#case Md287M5
if (((count(`i'5,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md156c
if (((count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))) Discard;
#break
#case Md207
#case MnonoM2
#case MnonoM5
#case MnonoM8
#case MnonoM11
#case MlalaM2
#case MlalaM5
#case Md308M5
#case Md378M6
#case Md378M11
#case Md354M8
#case Md385M3
#case Md385M8
#case Md310M7
#case Md315M5
#case Md377M3
#case Md377M8
if (((count(`i'8,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md207c
if (((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md378
#case Mno6M2
#case Mno6M6
#case Md388M3
#case Md387M3
#case Md397M5
#case Md401M3
#case Md399M4
#case Md399M5
if (((count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md270
#case Md308M8
#case Md354M3
#case Md310M8
#case Md323M5
#case Md307M3
#case Md307M4
#case Md307M7
#case Md307M8
#case Md313M8
if (((count(`i'8,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md270c
if (((count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md310
#case Md389M3
#case Md387M4
#case Md406M3
#case Md404M4
#case Md424M8
#case Md424M9
#case Md423M9
if (((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md310c
if (((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md218
#case MnonoM7
#case Md325M7
#case Md370M7
#case Md371M5
#case Md329M6
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md218c
if (((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md175
#case Md271d
#case Md195M5
#case Md195M6
if (((count(`i'5,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md175c
if (((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))) Discard;
#break
#case Md271
#case Md308M11
#case Md353M7
#case Md325M8
#case Md320M5
#case Md309M6
#case Md309M7
#case Md331M6
if (((count(`i'8,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md271c
if (((count(`i'5,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md325
#case Md390M3
#case Md390M8
#case Md387M5
#case Md407M5
#case Md406M5
#case Md420M4
#case Md420M5
#case Md420M7
#case Md420M8
#case Md419M6
if (((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md325c
if (((count(`i'6,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'2,1) >= 0))) Discard;
#break
#case Mfastar2
#case Mnostar5M2
#case Mlastar5M2
#case Md208d
#case Md349d
#case Md198M5
if (((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Mnostar5
#case Md389d
#case Md385d
if (((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md208
#case Md308M6
#case Md385M4
#case Md349M7
#case Md347M5
#case Md373M8
if (((count(`i'8,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md208c
if (((count(`i'4,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md385
#case Md389M8
#case Md387M6
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md385c
if (((count(`i'8,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))) Discard;
#break
#case Mnono
#case Md387M8
#case Md399M8
if (((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md308
#case Mno6M3
#case Mno6M7
#case Mno6M9
#case Mno6M12
#case Md388M8
#case Md387M9
#case Md396M5
#case Md406M9
#case Md424M4
#case Md424M5
if (((count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md387
#case Mj4M3
#case Mj4M9
#case Mno1M3
#case Mno1M7
#case Mno2M2
#case Mno2M3
#case Mno2M6
#case Mno2M7
if (((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))) Discard;
#break
#case Md138
#case Md286d
#case Md195M7
#case Md335d
#case Md319d
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md138c
if (((count(`i'3,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md286
#case Md309M8
#case Md335M6
#case Md307M6
#case Md319M8
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md286c
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md307
#case Md424M3
#case Md421M9
if (((count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md195
#case Md309d
#case Md418d
#case Md422d
if (((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md309
#case Md420M3
#case Md420M9
#case Md418M6
#case Md424M7
#case Md422M9
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md309c
if (((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md424
#case Mj1M10
#case MhahaM2
#case MhahaM5
#case Mno2M9
#case Mno2M13
if (((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md420
#case Mj6M7
#case MhahaM13
#case Mno2M10
if (((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))) Discard;
#break
#case Mno2
if (((count(`i'3,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))) Discard;
#break
#case Md173
#case Md234M5
#case Md256M1
#case Md256M7
#case Md266M1
#case Md264M6
#case Md275M5
#case Md255M6
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md173c
if (((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md90
#case Md115M1
#case Md114M3
#case Md174d
#case Md228d
if (((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md90c
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md174
#case Md219M1
#case Md232M2
#case Md256M3
#case Md273M5
#case Md252M5
#case Md231M7
#case Md226M1
#case Md228M3
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md174c
if (((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))) Discard;
#break
#case Md256
#case Md344M2
#case Md339M1
#case Md339M6
#case Md332M6
#case Md326M6
#case Md362M4
#case Md341M6
if (((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md256c
if (((count(`i'6,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md179
#case Md232M1
#case Md266M3
#case Md262M5
#case Md231M2
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md179c
if (((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'5,1) >= 0))) Discard;
#break
#case Md178
#case Md234M1
#case Md220M2
#case Md219M3
#case Md266M7
#case Md260M5
#case Md254M6
#case Md258M6
#case Md253M6
#case Md265M1
#case Md265M7
#case Md268M5
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md178c
if (((count(`i'5,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))) Discard;
#break
#case Md266
#case Md344M1
#case Md345M1
#case Md339M3
#case Md339M4
#case Md334M6
#case Md333M6
#case Md357M2
#case Md342M6
#case Md342M11
#case Md340M1
#case Md340M4
#case Md340M6
#case Md340M8
#case Md346M6
if (((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md266c
if (((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md69
#case Mfastar3M4
#case Md127M2
#case Md141M5
#case Md267d
#case Md125M6
#case Md253d
#case Md126M1
#case Md126M6
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md69c
if (((count(`i'4,1) >= 0))||((count(`i'3,1) >= 0))) Discard;
#break
#case Md64
#case Md127M3
#case Md273d
#case Md143M5
#case Md231d
#case Md131M1
if (((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md64c
if (((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Mt1star24
#case Md76d
#case Md142d
#case Md122d
if (((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md76
#case Mfastar3M3
#case Md127M6
#case Md187d
#case Md142M5
#case Md248d
#case Md122M6
if (((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md76c
if (((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md127
#case Mnostar6M3
#case Mnostar6M4
#case Md219d
#case Md327d
#case Md196M6
#case Md357d
#case Md199M1
if (((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md187
#case Md220M3
#case Md219M7
#case Md250M5
#case Md248M6
if (((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md187c
if (((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md219
#case Md379M1
#case Md379M2
#case Md345M2
#case Md339M5
#case Md339M8
#case Md327M6
#case Md328M6
#case Md357M7
#case Md343M1
if (((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md219c
if (((count(`i'7,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md339
#case Md391M1
#case Md391M2
#case Md386M7
#case Md416M1
if (((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md106
#case Md148M3
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md106c
if (((count(`i'3,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'1,1) >= 0))) Discard;
#break
#case Md55
#case Md149M3
#case Md149M5
#case Md148M4
#case Md290d
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md55c
if (((count(`i'2,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))) Discard;
#break
#case Md107
#case Md148M5
if (((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md107c
if (((count(`i'3,1) >= 0))||((count(`i'2,1) >= 0))) Discard;
#break
#case Md148
#case Md192M1
#case Md192M4
#case Md192M6
#case Md192M10
#case Md298d
if (((count(`i'2,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md148c
if (((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'1,1) >= 0))) Discard;
#break
#case Md54
#case Md149M4
#case Md291d
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md54c
if (((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))) Discard;
#break
#case Md149
#case Md192M3
#case Md192M5
#case Md294d
#case Md297d
if (((count(`i'1,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md149c
if (((count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))) Discard;
#break
#case Md192
#case Md340d
#case Md411d
#case Md426d
if (((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md67
#case Md177d
#case Md126M3
#case Md131M2
#case Md269d
if (((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md67c
if (((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'3,1) >= 0))) Discard;
#break
#case Md126
#case Md265d
#case Md358d
#case Md199M2
#case Md350d
#case Md355d
if (((count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md126c
if (((count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))) Discard;
#break
#case Md177
#case Md226M2
#case Md265M3
#case Md230M3
#case Md269M5
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md177c
if (((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md265
#case Md343M2
#case Md358M4
#case Md341M11
#case Md340M3
#case Md340M5
#case Md350M6
#case Md355M6
if (((count(`i'8,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md265c
if (((count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md340
#case Md416M2
#case Md411M7
if (((count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md340c
if (((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md186
#case Md235M6
#case Md234M10
#case Md232M7
#case Md244M4
#case Md244M6
#case Md246M6
#case Md233M3
#case Md233M7
#case Md247M6
#case Md251M6
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md186c
if (((count(`i'5,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))) Discard;
#break
#case Md234
#case MlalaM3
#case MlalaM6
#case Md344M8
#case Md345M6
#case Md354M7
#case Md299M7
#case Md303M7
#case Md342M5
#case Md342M8
#case Md341M1
#case Md341M7
#case Md348M6
if (((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md234c
if (((count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md88
#case Md114M2
#case Md230d
if (((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md88c
if (((count(`i'3,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md91
#case Md119M3
#case Md114M9
#case Md227d
if (((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0))) Discard;
#break
#case Md91c
if (((count(`i'4,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'1,1) >= 0))) Discard;
#break
#case Md114
#case Mlastar5M1
#case Mlastar5M3
#case Md226d
#case Md360d
#case Md376d
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md114c
if (((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md180
#case Md235M1
#case Md220M6
#case Md244M5
#case Md243M6
#case Md241M6
#case Md261M6
#case Md226M7
#case Md227M3
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md180c
if (((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md226
#case Md343M4
#case Md360M4
#case Md373M2
#case Md373M11
#case Md341M3
#case Md376M3
#case Md351M6
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md226c
if (((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md244
#case Md345M8
#case Md356M7
#case Md366M3
#case Md366M7
#case Md370M6
#case Md370M11
#case Md367M4
#case Md341M5
if (((count(`i'8,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md244c
if (((count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md341
#case Md406M2
#case Md403M4
#case Md416M4
#case Md412M7
if (((count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md342
#case Md399M1
#case Md398M2
#case Md416M5
#case Md413M7
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md342c
if (((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md70
#case Mfastar3M9
#case Md117M6
#case Md261d
#case Md131M6
if (((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md70c
if (((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md131
#case Md199M4
#case Md351d
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md131c
if (((count(`i'5,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Mfastar3
#case Mnostar6M2
#case Mnostar6M7
#case Md220d
#case Md353d
#case Md193M7
#case Md300d
#case Md199M7
if (((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md199
#case Md410d
#case Md343d
#case Md414d
if (((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md220
#case MlalaM1
#case MlalaM4
#case Md379M4
#case Md379M8
#case Md345M5
#case Md353M6
#case Md302M7
#case Md300M7
#case Md343M8
if (((count(`i'8,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md220c
if (((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'2,1) >= 0))) Discard;
#break
#case Md343
#case Md410M5
#case Md406M12
#case Md416M6
#case Md414M7
if (((count(`i'7,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md343c
if (((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md113
#case Md115M5
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md113c
if (((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md111
#case Md119M4
#case Md115M6
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md111c
if (((count(`i'4,1) >= 0))||((count(`i'2,1) >= 0))) Discard;
#break
#case Md115
#case Mlastar5M4
#case Mlastar5M6
#case Md232d
#case Md375d
#case Md364d
if (((count(`i'7,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md115c
if (((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md232
#case Md344M4
#case Md345M4
#case Md349M6
#case Md375M3
#case Md375M8
#case Md364M4
#case Md373M1
#case Md373M6
if (((count(`i'8,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md232c
if (((count(`i'5,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md344
#case Md391M5
#case Md391M7
#case Md389M7
#case Md409M5
#case Md406M7
#case Md416M7
#case Md415M7
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md344c
if (((count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md345
#case Md391M4
#case Md391M9
#case Md388M7
#case Md402M3
#case Md402M8
#case Md399M7
#case Md399M12
#case Md406M1
#case Md416M9
if (((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md416
#case Mj2M8
#case MhahaM1
#case MhahaM3
#case MhahaM4
#case MhahaM6
#case Mno1M1
#case Mno1M5
if (((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Mnostar6
#case Md379d
#case Md390d
#case Md398d
if (((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md379
#case Md391M6
#case Md390M7
#case Md399M2
#case Md398M8
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md379c
if (((count(`i'8,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))) Discard;
#break
#case Md189
#case Md235M10
#case Md245M6
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md189c
if (((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md235
#case MlalaM7
#case MlalaM11
#case Md295M7
#case Md296M7
#case Md370M5
#case Md370M8
#case Md369M4
#case Md373M5
#case Md373M7
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md235c
if (((count(`i'3,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))) Discard;
#break
#case Md370
#case Md400M3
#case Md400M8
#case Md399M3
#case Md408M5
#case Md406M6
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md370c
if (((count(`i'6,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Mlala
#case Md401M8
#case Md399M6
#case Md399M9
#case Md406M8
if (((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md399
#case Mj11M3
#case Mj11M9
#case Mno1M2
#case Mno1M6
if (((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md391
#case Mj4M8
#case Mno1M4
#case Mno1M8
if (((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md112
#case Md119M9
if (((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md112c
if (((count(`i'3,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'1,1) >= 0))) Discard;
#break
#case Md119
#case Mlastar5M7
#case Mlastar5M10
#case Md295d
#case Md369d
if (((count(`i'1,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md119c
if (((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'1,1) >= 0))) Discard;
#break
#case Mlastar5
#case Md373d
#case Md405d
if (((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md373
#case Md406M4
#case Md405M4
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md373c
if (((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md406
#case Mj5M5
#case MhahaM7
#case MhahaM8
#case MhahaM9
#case MhahaM10
#case Mno1M9
#case Mno1M10
if (((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))) Discard;
#break
#case Mno6
#case Mno3M5
#case Mno1M13
if (((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))) Discard;
#break
#case Mno1
if (((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))) Discard;
#break
#case Mhaha
if (((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md84
#case Md116M1
#case Md122M5
#case Md123M9
#case Md238d
if (((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))) Discard;
#break
#case Md84c
if (count(`i'4,1) >= 0) Discard;
#break
#case Md56
#case Md150M1
#case Md150M2
#case Md125M9
#case Md122M9
#case Md288d
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md56c
if (((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md122
#case Md193M1
#case Md365d
#case Md304d
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md122c
if (((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md82
#case Md116M2
#case Md125M5
#case Md123M5
#case Md240d
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md82c
if (((count(`i'4,1) >= 0))||((count(`i'3,1) >= 0))) Discard;
#break
#case Md125
#case Md193M2
#case Md361d
#case Md301d
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md125c
if (((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md85
#case Md116M4
#case Md117M5
#case Md120M5
#case Md120M9
#case Md237d
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md85c
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md105
#case Md150M4
#case Md117M9
if (((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0))) Discard;
#break
#case Md105c
if (((count(`i'3,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md117
#case Md193M4
#case Md371d
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md117c
if (((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'3,1) >= 0))) Discard;
#break
#case Md116
#case Md225d
#case Md193M6
#case Md299d
#case Md201M6
#case Md201M10
#case Md372d
if (((count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md116c
if (((count(`i'5,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))) Discard;
#break
#case Md150
#case Md289d
#case Md193M10
#case Md293d
#case Md292d
if (((count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md150c
if (((count(`i'5,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))) Discard;
#break
#case Md193
#case Md302d
#case Md407d
#case Md425d
if (((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md182
#case Md225M2
#case Md258M5
#case Md253M7
#case Md223M6
#case Md240M1
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md182c
if (((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md253
#case Md299M1
#case Md357M4
#case Md300M6
#case Md358M1
#case Md358M8
#case Md363M1
#case Md301M2
#case Md301M4
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))) Discard;
#break
#case Md253c
if (((count(`i'2,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md183
#case Md225M1
#case Md248M5
#case Md255M7
#case Md223M10
#case Md238M1
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md183c
if (((count(`i'5,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'2,1) >= 0))) Discard;
#break
#case Md255
#case Md299M2
#case Md357M6
#case Md303M6
#case Md359M1
#case Md362M1
#case Md362M7
#case Md304M2
#case Md304M11
if (((count(`i'1,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))) Discard;
#break
#case Md255c
if (((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'2,1) >= 0))) Discard;
#break
#case Md188
#case Md225M4
#case Md241M5
#case Md247M7
#case Md224M6
#case Md224M10
#case Md237M1
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md188c
if (((count(`i'3,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))) Discard;
#break
#case Md247
#case Md299M4
#case Md366M4
#case Md366M6
#case Md296M6
#case Md368M1
#case Md367M1
#case Md367M7
#case Md298M4
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))) Discard;
#break
#case Md247c
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md225
#case Md302M6
#case Md299M8
#case Md377M7
#case Md377M11
#case Md372M1
if (((count(`i'8,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md225c
if (((count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md299
#case Md402M5
#case Md402M7
#case Md401M7
#case Md404M1
#case Md403M1
#case Md403M8
#case Md425M2
if (((count(`i'9,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md299c
if (((count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md261
#case Md295M1
#case Md300M4
#case Md360M8
#case Md298M3
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md261c
if (((count(`i'2,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md152
#case Md289M1
#case Md289M2
#case Md258M7
#case Md248M7
#case Md249M10
#case Md288M6
#case Md288M7
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md152c
if (((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'1,1) >= 0))) Discard;
#break
#case Md249
#case Md293M1
#case Md293M2
#case Md303M5
#case Md300M5
#case Md361M8
#case Md365M7
#case Md292M4
#case Md292M8
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md249c
if (((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'5,1) >= 0))) Discard;
#break
#case Md223
#case Md377M1
#case Md377M2
#case Md303M8
#case Md300M8
#case Md359M8
#case Md363M7
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md223c
if (((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0))) Discard;
#break
#case Md248
#case Md302M1
#case Md357M5
#case Md300M11
#case Md365M1
#case Md304M4
#case Md304M8
if (((count(`i'8,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md248c
if (((count(`i'2,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md300
#case Md401M1
#case Md398M4
#case Md398M5
#case Md410M9
#case Md425M4
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md300c
if (((count(`i'6,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md151
#case Md289M4
#case Md241M7
#case Md242M10
#case Md290M7
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md151c
if (((count(`i'1,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))) Discard;
#break
#case Md290
#case Md294M6
#case Md294M11
#case Md292M5
#case Md298M5
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md290c
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md288
#case Md292M6
#case Md292M11
#case Md304M5
#case Md301M5
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))) Discard;
#break
#case Md288c
if (((count(`i'3,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0))) Discard;
#break
#case Md289
#case Md302M8
#case Md293M11
#case Md292M7
if (((count(`i'8,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md289c
if (((count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md292
#case Md427M1
#case Md427M2
#case Md427M7
#case Md427M12
#case Md425M5
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md292c
if (((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'3,1) >= 0))) Discard;
#break
#case Md251
#case Md295M2
#case Md303M4
#case Md364M7
#case Md298M7
if (((count(`i'1,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md251c
if (((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'5,1) >= 0))) Discard;
#break
#case Md241
#case Md302M4
#case Md366M5
#case Md366M8
#case Md296M11
#case Md371M1
#case Md298M8
if (((count(`i'8,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md241c
if (((count(`i'2,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md298
#case Md426M7
#case Md426M12
#case Md425M6
if (((count(`i'9,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md298c
if (((count(`i'4,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))) Discard;
#break
#case Md258
#case Md302M2
#case Md357M8
#case Md303M11
#case Md361M1
#case Md301M7
#case Md301M8
if (((count(`i'8,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md258c
if (((count(`i'2,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md301
#case Md425M7
if (((count(`i'9,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md301c
if (((count(`i'4,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))) Discard;
#break
#case Md302
#case Md402M6
#case Md402M9
#case Md401M12
#case Md407M1
#case Md425M8
if (((count(`i'9,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md302c
if (((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md303
#case Md401M2
#case Md398M6
#case Md398M9
#case Md409M8
#case Md425M9
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md303c
if (((count(`i'6,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md304
#case Md425M12
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md304c
if (((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))) Discard;
#break
#case Md425
#case Mla4M1
#case Mla4M2
#case Mla4M8
#case Mla4M13
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md425c
if (((count(`i'8,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))) Discard;
#break
#case Md245
#case Md295M4
#case Md296M4
#case Md369M7
#case Md297M3
#case Md297M4
#case Md297M7
#case Md297M8
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md245c
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0))) Discard;
#break
#case Md242
#case Md293M4
#case Md296M5
#case Md371M7
#case Md294M4
#case Md294M8
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md242c
if (((count(`i'3,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0))) Discard;
#break
#case Md224
#case Md377M4
#case Md296M8
#case Md368M7
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md224c
if (((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md296
#case Md401M4
#case Md400M4
#case Md400M5
#case Md400M6
#case Md400M9
#case Md408M8
#case Md426M4
#case Md426M9
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md296c
if (((count(`i'6,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md295
#case Md401M5
#case Md405M8
#case Md426M3
#case Md426M8
if (((count(`i'8,1) >= 0))||((count(`i'9,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md295c
if (((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))) Discard;
#break
#case Md293
#case Md401M6
#case Md407M8
#case Md427M3
#case Md427M5
#case Md427M8
#case Md427M9
if (((count(`i'8,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md293c
if (((count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md123
#case Md201M1
#case Md201M2
#case Md359d
#case Md363d
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md123c
if (((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md120
#case Md201M4
#case Md368d
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md120c
if (((count(`i'3,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md201
#case Md377d
#case Md404d
if (((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md377
#case Md401M9
#case Md404M8
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md377c
if (((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md401
#case Mj11M5
#case Mj11M6
#case Mj11M7
#case Mj11M10
#case Mj5M9
#case Mla4M3
#case Mla4M5
#case Mla4M9
#case Mla4M10
if (((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))) Discard;
#break
#case Md291
#case Md294M5
#case Md297M5
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md291c
if (((count(`i'5,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md294
#case Md427M4
#case Md427M6
#case Md426M5
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md294c
if (((count(`i'1,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md297
#case Md426M6
if (((count(`i'9,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md297c
if (((count(`i'4,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md426
#case Mla4M4
#case Mla4M7
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md426c
if (((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md427
#case Mla4M6
if (((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Mla4
if (((count(`i'7,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md80
#case Md121M2
#case Md191M1
#case Md243d
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md80c
if (((count(`i'2,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md78
#case Md121M3
#case Md121M4
#case Md191M9
#case Md190M1
#case Md190M9
#case Md246d
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md78c
if (((count(`i'4,1) >= 0))||((count(`i'1,1) >= 0))) Discard;
#break
#case Md121
#case Md321d
#case Md197M1
#case Md197M10
#case Md356d
#case Md367d
#case Md322d
if (((count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md121c
if (((count(`i'5,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'1,1) >= 0))) Discard;
#break
#case Md190
#case Md197M3
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md190c
if (((count(`i'3,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md191
#case Md197M4
#case Md197M5
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md191c
if (((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md197
#case Md336d
#case Md408d
#case Md417d
if (((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md108
#case Md124M3
#case Md143M9
if (((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))) Discard;
#break
#case Md108c
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md62
#case Md162d
#case Md124M4
#case Md124M5
#case Md142M9
#case Md141M9
#case Md275d
if (((count(`i'4,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))) Discard;
#break
#case Md62c
if (count(`i'2,1) >= 0) Discard;
#break
#case Md124
#case Md257d
#case Md196M10
#case Md326d
#case Md362d
if (((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md124c
if (((count(`i'5,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'1,1) >= 0))) Discard;
#break
#case Md158
#case Md285M1
#case Md229M2
#case Md257M3
#case Md273M10
#case Md252M10
#case Md228M7
#case Md274M1
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md158c
if (((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'5,1) >= 0))) Discard;
#break
#case Md162
#case Md239M2
#case Md285M10
#case Md284M1
#case Md284M10
#case Md257M4
#case Md257M5
#case Md250M10
#case Md254M10
#case Md275M10
#case Md238M7
#case Md277M1
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md162c
if (((count(`i'4,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))) Discard;
#break
#case Md257
#case Md321M2
#case Md336M1
#case Md336M11
#case Md328M11
#case Md326M11
#case Md362M8
#case Md322M1
if (((count(`i'8,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md257c
if (((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md284
#case Md311M1
#case Md336M4
#case Md334M11
#case Md317M1
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md284c
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'2,1) >= 0))) Discard;
#break
#case Md285
#case Md316M1
#case Md316M2
#case Md311M2
#case Md336M5
#case Md336M6
#case Md327M11
#case Md332M11
#case Md312M1
#case Md318M1
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md285c
if (((count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md336
#case Md392M1
#case Md392M2
#case Md386M12
#case Md417M1
if (((count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md336c
if (((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md229
#case Md311M4
#case Md321M4
#case Md349M11
#case Md376M8
#case Md312M11
#case Md319M1
#case Md319M5
if (((count(`i'1,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md229c
if (((count(`i'3,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md239
#case Md316M4
#case Md316M11
#case Md311M11
#case Md321M5
#case Md321M6
#case Md353M11
#case Md354M11
#case Md372M8
#case Md318M11
#case Md317M11
#case Md313M1
#case Md313M5
if (((count(`i'1,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))) Discard;
#break
#case Md239c
if (((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0))) Discard;
#break
#case Md321
#case Md392M4
#case Md392M12
#case Md388M12
#case Md403M9
#case Md417M12
#case Md421M1
#case Md421M5
if (((count(`i'9,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md321c
if (((count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md311
#case Md392M5
#case Md392M7
#case Md389M12
#case Md423M1
#case Md422M1
#case Md422M6
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))) Discard;
#break
#case Md311c
if (((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))) Discard;
#break
#case Md316
#case Md392M6
#case Md390M12
#case Md423M6
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md316c
if (((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md392
#case Mj4M13
#case Mj1M1
#case Mj1M6
if (((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md283
#case Md320M3
#case Md337M5
#case Md318M3
if (((count(`i'1,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md283c
if (((count(`i'6,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'1,1) >= 0))) Discard;
#break
#case Md61
#case Md154d
#case Md145M1
#case Md145M3
#case Md133M1
#case Md134M9
#case Md147M6
#case Md147M9
#case Md279d
#case Md277d
if (((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))) Discard;
#break
#case Md61c
if (((count(`i'2,1) >= 0))||((count(`i'1,1) >= 0))) Discard;
#break
#case Md154
#case Md287M1
#case Md287M3
#case Md276M1
#case Md278M7
#case Md282M6
#case Md282M7
#case Md279M7
#case Md277M3
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md154c
if (((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))) Discard;
#break
#case Md279
#case Md314M3
#case Md338M7
#case Md323M1
#case Md323M3
#case Md337M6
#case Md337M8
#case Md317M3
#case Md322M3
#case Md322M8
#case Md346M7
if (((count(`i'1,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md279c
if (((count(`i'5,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'1,1) >= 0))) Discard;
#break
#case Md109
#case Md133M3
#case Md134M6
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md109c
if (((count(`i'4,1) >= 0))||((count(`i'1,1) >= 0))) Discard;
#break
#case Md134
#case Md278d
#case Md194M3
#case Md198M1
#case Md198M3
#case Md312d
#case Md346d
if (((count(`i'1,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md134c
if (((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md169
#case Md276M3
#case Md278M6
#case Md267M7
#case Md268M2
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md169c
if (((count(`i'3,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md278
#case Md315M3
#case Md338M8
#case Md347M1
#case Md347M3
#case Md337M7
#case Md337M11
#case Md312M3
#case Md346M2
if (((count(`i'8,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md278c
if (((count(`i'1,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md337
#case Md396M3
#case Md393M7
#case Md393M8
#case Md417M3
if (((count(`i'1,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md337c
if (((count(`i'5,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md277
#case Md312M4
#case Md317M4
#case Md322M4
#case Md322M6
#case Md348M11
#case Md313M2
#case Md313M11
if (((count(`i'2,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md277c
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'2,1) >= 0))) Discard;
#break
#case Md274
#case Md318M4
#case Md322M5
#case Md351M11
#case Md319M2
#case Md319M11
if (((count(`i'2,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md274c
if (((count(`i'6,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md322
#case Md417M4
#case Md412M12
#case Md421M2
#case Md421M12
if (((count(`i'9,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md322c
if (((count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md147
#case Md282d
#case Md194M1
#case Md317d
if (((count(`i'2,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md147c
if (((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))) Discard;
#break
#case Md282
#case Md315M1
#case Md338M6
#case Md317M8
if (((count(`i'8,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md282c
if (((count(`i'2,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md317
#case Md417M5
#case Md413M12
#case Md423M12
if (((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md317c
if (((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md272
#case Md352M1
#case Md352M3
#case Md338M11
#case Md320M1
#case Md318M8
#case Md350M2
#case Md350M7
if (((count(`i'2,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md272c
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md318
#case Md417M6
#case Md414M12
#case Md422M2
#case Md422M12
if (((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md318c
if (((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md267
#case Md314M1
#case Md338M5
#case Md312M8
#case Md355M2
#case Md355M7
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md267c
if (((count(`i'3,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md312
#case Md417M7
#case Md415M12
#case Md423M2
if (((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md312c
if (((count(`i'2,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md338
#case Md397M1
#case Md397M3
#case Md396M1
#case Md393M5
#case Md393M6
#case Md393M9
#case Md393M12
#case Md417M9
#case Md411M2
#case Md411M8
if (((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md338c
if (((count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0))) Discard;
#break
#case Md417
#case Mj2M13
#case Mj1M2
#case Mj1M13
if (((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md417c
if (((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md145
#case Md287d
#case Md194M10
#case Md323d
#case Md313d
if (((count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md145c
if (((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))) Discard;
#break
#case Md287
#case Md315M8
#case Md323M11
#case Md313M4
if (((count(`i'7,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md287c
if (((count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md323
#case Md396M7
#case Md413M3
#case Md421M3
#case Md421M4
if (((count(`i'7,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md323c
if (((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md319
#case Md422M3
#case Md421M6
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md319c
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md313
#case Md423M3
#case Md421M7
#case Md421M8
if (((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md313c
if (((count(`i'3,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0))) Discard;
#break
#case Md421
#case Mj1M3
if (((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md320
#case Md396M6
#case Md414M3
#case Md422M4
#case Md422M5
if (((count(`i'8,1) >= 0))||((count(`i'9,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))) Discard;
#break
#case Md320c
if (((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))) Discard;
#break
#case Md133
#case Md276d
#case Md194M7
#case Md314d
#case Md198M7
#case Md198M10
#case Md348d
if (((count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md133c
if (((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md198
#case Md347d
#case Md415d
if (((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md276
#case Md315M7
#case Md314M8
#case Md347M8
#case Md347M11
#case Md348M3
if (((count(`i'7,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md276c
if (((count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md347
#case Md396M8
#case Md415M3
if (((count(`i'7,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md347c
if (((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md314
#case Md397M6
#case Md397M8
#case Md396M9
#case Md415M8
#case Md412M3
#case Md412M8
#case Md423M4
if (((count(`i'7,1) >= 0))||((count(`i'9,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md314c
if (((count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md194
#case Md315d
#case Md413d
#case Md423d
if (((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md315
#case Md397M7
#case Md397M9
#case Md396M12
#case Md413M8
#case Md423M5
if (((count(`i'9,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md315c
if (((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md396
#case Mno3M6
#case Mno3M7
#case Mno3M8
#case Mno3M9
#case Mj2M3
#case Mj1M4
#case Mj1M5
if (((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))) Discard;
#break
#case Md422
#case Mj1M7
if (((count(`i'7,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md422c
if (((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md423
#case Mj1M8
#case Mj1M9
if (((count(`i'7,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md423c
if (((count(`i'8,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Mj1
if (((count(`i'2,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))) Discard;
#break
#case Md275
#case Md332M1
#case Md334M1
#case Md326M3
#case Md326M5
#case Md348M1
#case Md331M1
#case Md330M1
#case Md330M11
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md275c
if (((count(`i'3,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))) Discard;
#break
#case Md273
#case Md327M1
#case Md326M4
#case Md351M1
#case Md329M1
#case Md329M2
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md273c
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'3,1) >= 0))) Discard;
#break
#case Md326
#case Md386M1
#case Md412M1
#case Md419M1
if (((count(`i'8,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md326c
if (((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md262
#case Md349M1
#case Md334M3
#case Md333M3
#case Md360M11
#case Md346M8
#case Md335M1
#case Md335M2
#case Md335M5
#case Md335M11
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md262c
if (((count(`i'1,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))) Discard;
#break
#case Md260
#case Md353M1
#case Md327M3
#case Md333M4
#case Md361M11
#case Md350M1
#case Md350M8
#case Md331M2
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md260c
if (((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))) Discard;
#break
#case Md264
#case Md354M2
#case Md332M8
#case Md333M5
#case Md359M11
#case Md346M1
#case Md331M11
if (((count(`i'1,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))) Discard;
#break
#case Md264c
if (((count(`i'3,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md333
#case Md388M1
#case Md386M3
#case Md410M12
#case Md411M1
#case Md411M9
#case Md418M1
#case Md418M5
if (((count(`i'8,1) >= 0))||((count(`i'9,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md333c
if (((count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))) Discard;
#break
#case Md141
#case Md254d
#case Md196M5
#case Md331d
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md141c
if (((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md254
#case Md354M1
#case Md334M8
#case Md328M5
#case Md363M8
#case Md331M5
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md254c
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md334
#case Md389M1
#case Md386M4
#case Md413M1
#case Md418M2
#case Md418M12
if (((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md334c
if (((count(`i'3,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md142
#case Md250d
#case Md196M4
#case Md330d
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md142c
if (((count(`i'5,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'2,1) >= 0))) Discard;
#break
#case Md250
#case Md353M2
#case Md327M8
#case Md328M4
#case Md365M8
#case Md330M2
#case Md330M5
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md250c
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md327
#case Md390M1
#case Md390M2
#case Md386M5
#case Md414M1
#case Md419M2
if (((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md327c
if (((count(`i'7,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md143
#case Md252d
#case Md196M3
#case Md329d
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md143c
if (((count(`i'5,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'2,1) >= 0))) Discard;
#break
#case Md252
#case Md349M2
#case Md332M3
#case Md328M3
#case Md364M8
#case Md329M5
#case Md329M11
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md252c
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md332
#case Md389M2
#case Md386M6
#case Md415M1
#case Md419M12
if (((count(`i'5,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md332c
if (((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'2,1) >= 0))) Discard;
#break
#case Md196
#case Md328d
#case Md409d
#case Md419d
if (((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md328
#case Md388M2
#case Md386M9
#case Md409M9
#case Md419M5
if (((count(`i'1,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md328c
if (((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md386
#case Mj4M1
#case Mj4M2
#case Mj2M1
#case Mj6M1
#case Mj6M2
#case Mj6M6
#case Mj6M13
if (((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))) Discard;
#break
#case Md331
#case Md419M9
#case Md418M7
#case Md418M8
if (((count(`i'8,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md331c
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'2,1) >= 0))) Discard;
#break
#case Md335
#case Md418M9
if (((count(`i'8,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md335c
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))) Discard;
#break
#case Md418
#case Mj6M3
#case Mj6M10
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md418c
if (((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md329
#case Md419M7
if (((count(`i'8,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))) Discard;
#break
#case Md329c
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'2,1) >= 0))) Discard;
#break
#case Md330
#case Md419M8
if (((count(`i'8,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md330c
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md419
#case Mj6M4
#case Mj6M5
#case Mj6M8
#case Mj6M9
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md419c
if (((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))) Discard;
#break
#case Mj6
if (((count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md268
#case Md348M2
#case Md346M3
#case Md355M3
#case Md355M5
if (((count(`i'1,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md268c
if (((count(`i'5,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md269
#case Md351M2
#case Md350M3
#case Md355M4
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md269c
if (((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'3,1) >= 0))) Discard;
#break
#case Md355
#case Md412M2
#case Md411M3
if (((count(`i'8,1) >= 0))||((count(`i'9,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md355c
if (((count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md346
#case Md415M2
#case Md413M2
#case Md411M4
#case Md411M6
if (((count(`i'9,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md346c
if (((count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'1,1) >= 0))) Discard;
#break
#case Md350
#case Md414M2
#case Md411M5
if (((count(`i'9,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md350c
if (((count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md411
#case Mj2M2
if (((count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md411c
if (((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md246
#case Md354M4
#case Md349M8
#case Md356M4
#case Md356M6
#case Md369M8
#case Md368M8
#case Md348M8
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))) Discard;
#break
#case Md246c
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md348
#case Md415M4
#case Md413M4
#case Md412M4
#case Md412M6
if (((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md348c
if (((count(`i'5,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md243
#case Md353M8
#case Md356M5
#case Md371M8
#case Md351M8
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md243c
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md351
#case Md414M4
#case Md412M5
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md351c
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md356
#case Md388M9
#case Md408M9
#case Md412M9
if (((count(`i'9,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md356c
if (((count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md412
#case Mj2M4
if (((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md354
#case Md389M9
#case Md388M6
#case Md404M9
#case Md413M9
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md354c
if (((count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'3,1) >= 0))) Discard;
#break
#case Md413
#case Mj2M5
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))) Discard;
#break
#case Md413c
if (((count(`i'8,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md352
#case Md397M12
#case Md414M8
if (((count(`i'8,1) >= 0))||((count(`i'9,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md352c
if (((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))) Discard;
#break
#case Md353
#case Md390M4
#case Md390M9
#case Md388M5
#case Md407M9
#case Md414M9
if (((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md353c
if (((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))) Discard;
#break
#case Md414
#case Mj2M6
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md414c
if (((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md349
#case Md389M4
#case Md388M4
#case Md405M9
#case Md415M9
if (((count(`i'5,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md349c
if (((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))) Discard;
#break
#case Md415
#case Mj2M7
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md415c
if (((count(`i'8,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md397
#case Mno3M10
#case Mno3M13
#case Mj2M9
if (((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md388
#case Mj4M4
#case Mj4M10
#case Mj5M10
#case Mj2M10
if (((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Mj2
if (((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))) Discard;
#break
#case Md231
#case Md375M1
#case Md375M2
#case Md357M11
#case Md360M1
#case Md364M1
if (((count(`i'1,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md231c
if (((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md357
#case Md402M1
#case Md402M2
#case Md398M7
#case Md398M12
#case Md410M1
#case Md409M1
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md357c
if (((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'3,1) >= 0))) Discard;
#break
#case Md233
#case Md375M4
#case Md366M11
#case Md369M1
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md233c
if (((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'3,1) >= 0))) Discard;
#break
#case Md366
#case Md402M4
#case Md400M7
#case Md400M12
#case Md408M1
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md366c
if (((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md375
#case Md402M12
#case Md405M1
if (((count(`i'9,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md375c
if (((count(`i'7,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md402
#case Mj11M8
#case Mj11M13
#case Mj5M1
if (((count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md237
#case Md371M2
#case Md368M2
#case Md367M2
#case Md367M5
#case Md372M5
if (((count(`i'2,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md237c
if (((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md238
#case Md365M2
#case Md363M2
#case Md362M2
#case Md362M5
#case Md372M6
if (((count(`i'2,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md238c
if (((count(`i'4,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md240
#case Md361M2
#case Md359M2
#case Md358M2
#case Md358M5
#case Md372M11
if (((count(`i'2,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md240c
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'5,1) >= 0))) Discard;
#break
#case Md372
#case Md407M2
#case Md404M2
#case Md403M2
#case Md403M5
if (((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md372c
if (((count(`i'5,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'3,1) >= 0))) Discard;
#break
#case Md227
#case Md369M2
#case Md367M3
#case Md376M5
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md227c
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md228
#case Md364M2
#case Md362M3
#case Md376M6
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md228c
if (((count(`i'4,1) >= 0))||((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0))) Discard;
#break
#case Md230
#case Md360M2
#case Md358M3
#case Md376M11
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))) Discard;
#break
#case Md230c
if (((count(`i'1,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'4,1) >= 0))) Discard;
#break
#case Md376
#case Md405M2
#case Md403M3
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md376c
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md367
#case Md408M2
#case Md403M6
if (((count(`i'8,1) >= 0))||((count(`i'9,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md367c
if (((count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md362
#case Md409M2
#case Md403M7
if (((count(`i'9,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md362c
if (((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md358
#case Md410M2
#case Md403M12
if (((count(`i'7,1) >= 0))||((count(`i'9,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md358c
if (((count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md403
#case Mj5M2
if (((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md368
#case Md408M3
#case Md404M6
if (((count(`i'9,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))) Discard;
#break
#case Md368c
if (((count(`i'4,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'1,1) >= 0))) Discard;
#break
#case Md363
#case Md409M3
#case Md404M7
if (((count(`i'9,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))) Discard;
#break
#case Md363c
if (((count(`i'4,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'1,1) >= 0))) Discard;
#break
#case Md359
#case Md410M3
#case Md404M12
if (((count(`i'3,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md359c
if (((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'2,1) >= 0))) Discard;
#break
#case Md404
#case Mj5M3
if (((count(`i'3,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md404c
if (((count(`i'5,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'2,1) >= 0))) Discard;
#break
#case Md369
#case Md408M4
#case Md405M6
if (((count(`i'9,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md369c
if (((count(`i'4,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))) Discard;
#break
#case Md364
#case Md409M4
#case Md405M7
if (((count(`i'9,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md364c
if (((count(`i'7,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md360
#case Md410M4
#case Md405M12
if (((count(`i'3,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md360c
if (((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md405
#case Mj5M4
if (((count(`i'3,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))) Discard;
#break
#case Md405c
if (((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))) Discard;
#break
#case Md371
#case Md408M6
#case Md407M6
if (((count(`i'9,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md371c
if (((count(`i'4,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md365
#case Md409M6
#case Md407M7
if (((count(`i'9,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md365c
if (((count(`i'4,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md361
#case Md410M6
#case Md407M12
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md361c
if (((count(`i'4,1) >= 0))||((count(`i'5,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md407
#case Mj5M6
if (((count(`i'3,1) >= 0))||((count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md407c
if (((count(`i'5,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md408
#case Mj5M7
if (((count(`i'7,1) >= 0))||((count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md408c
if (((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md409
#case Mj5M8
if (((count(`i'7,1) >= 0))||((count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))) Discard;
#break
#case Md409c
if (((count(`i'6,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))) Discard;
#break
#case Md410
#case Mj5M13
if (((count(`i'7,1) >= 0))||((count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md410c
if (((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Mj5
if (((count(`i'7,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))) Discard;
#break
#case Md398
#case Mj11M1
#case Mj11M2
if (((count(`i'1,1) >= 0))||((count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'10,1) >= 0))) Discard;
#break
#case Md398c
if (((count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))) Discard;
#break
#case Md400
#case Mj11M4
if (((count(`i'1,1) >= 0))||((count(`i'2,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))) Discard;
#break
#case Md400c
if (((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Mj11
if (((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))) Discard;
#break
#case Md393
#case Mno3M1
#case Mno3M3
if (((count(`i'2,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Md393c
if (((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))) Discard;
#break
#case Mno3
if (((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'11,1) >= 0))) Discard;
#break
#case Md389
#case Mj4M5
#case Mj4M7
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'10,1) >= 0))) Discard;
#break
#case Md389c
if (((count(`i'8,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))) Discard;
#break
#case Md390
#case Mj4M6
if (((count(`i'5,1) >= 0))||((count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))) Discard;
#break
#case Md390c
if (((count(`i'8,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'3,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0))) Discard;
#break
#case Mj4
if (((count(`i'6,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'7,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'7,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0))||((count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'7,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'2,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'8,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'4,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'1,1) >= 0)&&(count(`i'3,1) >= 0)&&(count(`i'5,1) >= 0)&&(count(`i'8,1) >= 0))||((count(`i'3,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'4,1) >= 0)&&(count(`i'6,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0))||((count(`i'8,1) >= 0)&&(count(`i'9,1) >= 0)&&(count(`i'10,1) >= 0)&&(count(`i'11,1) >= 0))) Discard;
#break
#endswitch
endif;
#endprocedure
