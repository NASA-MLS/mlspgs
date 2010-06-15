module cfm_lsf_m ! This module is to help cfm_mlssetup subroutine.
                 ! Quantities allocated in here should be deallocate by cfm_mlssetup_m
    use MLSCommon, only: r8

    implicit none

    private

    public :: CreateLimbSidebandFractions, FillLimbSidebandFractions

!---------------------------- RCS Ident Info -------------------------------
   character(len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
   private :: not_used_here
!---------------------------------------------------------------------------

    real(r8), dimension(25), parameter :: vlimbSidebandFraction2L = &
    (/0.54213_r8, 0.53927_r8, 0.53572_r8, 0.53244_r8, 0.52971_r8, 0.52693_r8, 0.52451_r8, &
      0.52279_r8, 0.52160_r8, 0.52076_r8, 0.52018_r8, 0.51977_r8, 0.51948_r8, 0.51919_r8, &
      0.51878_r8, 0.51822_r8, 0.51742_r8, 0.51632_r8, 0.51479_r8, 0.51273_r8, 0.51052_r8, &
      0.50846_r8, 0.50613_r8, 0.50368_r8, 0.50160_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction3L = &
    (/0.52287_r8, 0.51885_r8, 0.51507_r8, 0.51216_r8, 0.50998_r8, 0.50797_r8, 0.50635_r8, &
      0.50526_r8, 0.50455_r8, 0.50405_r8, 0.50372_r8, 0.50348_r8, 0.50332_r8, 0.50316_r8, &
      0.50293_r8, 0.50263_r8, 0.50220_r8, 0.50163_r8, 0.50088_r8, 0.49992_r8, 0.49898_r8, &
      0.49819_r8, 0.49741_r8, 0.49676_r8, 0.49637_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction7L = &
    (/0.48358_r8, 0.48470_r8, 0.48580_r8, 0.48669_r8, 0.48739_r8, 0.48807_r8, 0.48865_r8, &
      0.48905_r8, 0.48932_r8, 0.48950_r8, 0.48964_r8, 0.48973_r8, 0.48979_r8, 0.48985_r8, &
      0.48995_r8, 0.49007_r8, 0.49024_r8, 0.49047_r8, 0.49078_r8, 0.49116_r8, 0.49153_r8, &
      0.49178_r8, 0.49191_r8, 0.49171_r8, 0.49100_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction8L = &
    (/0.47597_r8, 0.47481_r8, 0.47369_r8, 0.47277_r8, &
      0.47205_r8, 0.47135_r8, 0.47076_r8, 0.47034_r8, 0.47006_r8, &
      0.46987_r8, 0.46972_r8, 0.46962_r8, 0.46956_r8, 0.46949_r8, &
      0.46940_r8, 0.46927_r8, 0.46907_r8, 0.46882_r8, 0.46846_r8, &
      0.46799_r8, 0.46747_r8, 0.46700_r8, 0.46645_r8, 0.46589_r8, &
      0.46543_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction9L = &
    (/0.47767_r8, 0.47701_r8, 0.47627_r8, 0.47565_r8, &
      0.47517_r8, 0.47471_r8, 0.47433_r8, 0.47408_r8, 0.47390_r8, &
      0.47379_r8, 0.47371_r8, 0.47365_r8, 0.47361_r8, 0.47357_r8, &
      0.47352_r8, 0.47344_r8, 0.47333_r8, 0.47319_r8, 0.47299_r8, &
      0.47274_r8, 0.47248_r8, 0.47224_r8, 0.47195_r8, 0.47165_r8, &
      0.47136_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction15L = &
    (/0.95750_r8, 0.95800_r8, 0.95850_r8, 0.95890_r8, &
      0.95920_r8, 0.95960_r8, 0.95990_r8, 0.96010_r8, 0.96020_r8, &
      0.96030_r8, 0.96040_r8, 0.96040_r8, 0.96050_r8, 0.96050_r8, &
      0.96060_r8, 0.96060_r8, 0.96070_r8, 0.96090_r8, 0.96110_r8, &
      0.96140_r8, 0.96170_r8, 0.96210_r8, 0.96250_r8, 0.96300_r8, &
      0.96350_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction16L = &
    (/0.80400_r8, 0.80870_r8, 0.81290_r8, 0.81670_r8, &
      0.81960_r8, 0.82270_r8, 0.82520_r8, 0.82710_r8, 0.82840_r8, &
      0.82940_r8, 0.83000_r8, 0.83050_r8, 0.83080_r8, 0.83110_r8, &
      0.83160_r8, 0.83220_r8, 0.83310_r8, 0.83440_r8, 0.83630_r8, &
      0.83890_r8, 0.84170_r8, 0.84490_r8, 0.84840_r8, 0.85320_r8, &
      0.85720_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction17L = &
    (/0.42500_r8, 0.41920_r8, 0.41350_r8, 0.40850_r8, &
      0.40450_r8, 0.40030_r8, 0.39680_r8, 0.39440_r8, 0.39270_r8, &
      0.39150_r8, 0.39060_r8, 0.39000_r8, 0.38950_r8, 0.38910_r8, &
      0.38850_r8, 0.38760_r8, 0.38640_r8, 0.38460_r8, 0.38220_r8, &
      0.37860_r8, 0.37490_r8, 0.37090_r8, 0.36590_r8, 0.36010_r8, &
      0.35400_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction18L = &
    (/0.94650_r8, 0.94680_r8, 0.94700_r8, 0.94730_r8, &
      0.94740_r8, 0.94760_r8, 0.94780_r8, 0.94790_r8, 0.94800_r8, &
      0.94800_r8, 0.94800_r8, 0.94810_r8, 0.94810_r8, 0.94810_r8, &
      0.94810_r8, 0.94820_r8, 0.94820_r8, 0.94830_r8, 0.94840_r8, &
      0.94860_r8, 0.94870_r8, 0.94890_r8, 0.94910_r8, 0.94940_r8, &
      0.94970_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction19L = &
    (/0.78500_r8, 0.78980_r8, 0.79410_r8, 0.79800_r8, &
      0.80110_r8, 0.80400_r8, 0.80670_r8, 0.80860_r8, 0.80990_r8, &
      0.81090_r8, 0.81160_r8, 0.81200_r8, 0.81240_r8, 0.81270_r8, &
      0.81320_r8, 0.81380_r8, 0.81470_r8, 0.81600_r8, 0.81810_r8, &
      0.82060_r8, 0.82370_r8, 0.82670_r8, 0.83040_r8, 0.83510_r8, &
      0.83960_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction20L = &
    (/0.41870_r8, 0.41310_r8, 0.40650_r8, 0.40200_r8, &
      0.39800_r8, 0.39370_r8, 0.39040_r8, 0.38790_r8, 0.38620_r8, &
      0.38500_r8, 0.38410_r8, 0.38350_r8, 0.38300_r8, 0.38260_r8, &
      0.38200_r8, 0.38120_r8, 0.37990_r8, 0.37810_r8, 0.37570_r8, &
      0.37200_r8, 0.36850_r8, 0.36450_r8, 0.35930_r8, 0.35380_r8, &
      0.34760_r8/)
    real(r8), dimension(129), parameter :: vlimbSidebandFraction23L = &
    (/0.51922_r8, 0.51923_r8, 0.51923_r8, 0.51923_r8, &
      0.51924_r8, 0.51924_r8, 0.51924_r8, 0.51925_r8, 0.51925_r8, &
      0.51925_r8, 0.51926_r8, 0.51926_r8, 0.51926_r8, 0.51927_r8, &
      0.51927_r8, 0.51927_r8, 0.51928_r8, 0.51928_r8, 0.51928_r8, &
      0.51929_r8, 0.51929_r8, 0.51929_r8, 0.51930_r8, 0.51930_r8, &
      0.51930_r8, 0.51931_r8, 0.51931_r8, 0.51931_r8, 0.51932_r8, &
      0.51932_r8, 0.51932_r8, 0.51933_r8, 0.51933_r8, 0.51933_r8, &
      0.51934_r8, 0.51934_r8, 0.51934_r8, 0.51935_r8, 0.51935_r8, &
      0.51935_r8, 0.51936_r8, 0.51936_r8, 0.51936_r8, 0.51937_r8, &
      0.51937_r8, 0.51937_r8, 0.51938_r8, 0.51938_r8, 0.51938_r8, &
      0.51939_r8, 0.51939_r8, 0.51939_r8, 0.51940_r8, 0.51940_r8, &
      0.51940_r8, 0.51941_r8, 0.51941_r8, 0.51941_r8, 0.51942_r8, &
      0.51942_r8, 0.51942_r8, 0.51943_r8, 0.51943_r8, 0.51943_r8, &
      0.51944_r8, 0.51944_r8, 0.51944_r8, 0.51945_r8, 0.51945_r8, &
      0.51945_r8, 0.51946_r8, 0.51946_r8, 0.51946_r8, 0.51947_r8, &
      0.51947_r8, 0.51947_r8, 0.51948_r8, 0.51948_r8, 0.51948_r8, &
      0.51949_r8, 0.51949_r8, 0.51949_r8, 0.51950_r8, 0.51950_r8, &
      0.51950_r8, 0.51951_r8, 0.51951_r8, 0.51951_r8, 0.51952_r8, &
      0.51952_r8, 0.51952_r8, 0.51953_r8, 0.51953_r8, 0.51953_r8, &
      0.51954_r8, 0.51954_r8, 0.51954_r8, 0.51955_r8, 0.51955_r8, &
      0.51955_r8, 0.51956_r8, 0.51956_r8, 0.51957_r8, 0.51957_r8, &
      0.51957_r8, 0.51958_r8, 0.51958_r8, 0.51958_r8, 0.51959_r8, &
      0.51959_r8, 0.51959_r8, 0.51960_r8, 0.51960_r8, 0.51960_r8, &
      0.51961_r8, 0.51961_r8, 0.51961_r8, 0.51962_r8, 0.51962_r8, &
      0.51962_r8, 0.51963_r8, 0.51963_r8, 0.51963_r8, 0.51964_r8, &
      0.51964_r8, 0.51964_r8, 0.51965_r8, 0.51965_r8, 0.51965_r8/)
    real(r8), dimension(129), parameter :: vlimbSidebandFraction24L = &
    (/0.48982_r8, 0.48982_r8, 0.48982_r8, 0.48982_r8, &
      0.48982_r8, 0.48982_r8, 0.48982_r8, 0.48982_r8, 0.48982_r8, &
      0.48982_r8, 0.48982_r8, 0.48982_r8, 0.48982_r8, 0.48982_r8, &
      0.48982_r8, 0.48982_r8, 0.48982_r8, 0.48981_r8, 0.48981_r8, &
      0.48981_r8, 0.48981_r8, 0.48981_r8, 0.48981_r8, 0.48981_r8, &
      0.48981_r8, 0.48981_r8, 0.48981_r8, 0.48981_r8, 0.48981_r8, &
      0.48981_r8, 0.48981_r8, 0.48981_r8, 0.48981_r8, 0.48981_r8, &
      0.48981_r8, 0.48981_r8, 0.48981_r8, 0.48981_r8, 0.48981_r8, &
      0.48981_r8, 0.48980_r8, 0.48980_r8, 0.48980_r8, 0.48980_r8, &
      0.48980_r8, 0.48980_r8, 0.48980_r8, 0.48980_r8, 0.48980_r8, &
      0.48980_r8, 0.48980_r8, 0.48980_r8, 0.48980_r8, 0.48980_r8, &
      0.48980_r8, 0.48980_r8, 0.48980_r8, 0.48980_r8, 0.48980_r8, &
      0.48980_r8, 0.48980_r8, 0.48980_r8, 0.48980_r8, 0.48979_r8, &
      0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, &
      0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, &
      0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, &
      0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, &
      0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, &
      0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, &
      0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, &
      0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, &
      0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, 0.48979_r8, &
      0.48978_r8, 0.48978_r8, 0.48978_r8, 0.48978_r8, 0.48978_r8, &
      0.48978_r8, 0.48978_r8, 0.48978_r8, 0.48978_r8, 0.48978_r8, &
      0.48978_r8, 0.48978_r8, 0.48978_r8, 0.48978_r8, 0.48978_r8, &
      0.48978_r8, 0.48978_r8, 0.48978_r8, 0.48978_r8, 0.48978_r8/)
    real(r8), dimension(129), parameter :: vlimbSidebandFraction25L = &
    (/0.47355_r8, 0.47355_r8, 0.47355_r8, 0.47355_r8, 0.47355_r8, 0.47355_r8, 0.47355_r8, &
      0.47355_r8, 0.47355_r8, 0.47355_r8, 0.47356_r8, 0.47356_r8, 0.47356_r8, 0.47356_r8, &
      0.47356_r8, 0.47356_r8, 0.47356_r8, 0.47356_r8, 0.47356_r8, 0.47356_r8, 0.47356_r8, &
      0.47356_r8, 0.47356_r8, 0.47356_r8, 0.47357_r8, 0.47357_r8, 0.47357_r8, 0.47357_r8, &
      0.47357_r8, 0.47357_r8, 0.47357_r8, 0.47357_r8, 0.47357_r8, 0.47357_r8, 0.47357_r8, &
      0.47357_r8, 0.47357_r8, 0.47358_r8, 0.47358_r8, 0.47358_r8, 0.47358_r8, 0.47358_r8, &
      0.47358_r8, 0.47358_r8, 0.47358_r8, 0.47358_r8, 0.47358_r8, 0.47358_r8, 0.47358_r8, &
      0.47358_r8, 0.47359_r8, 0.47359_r8, 0.47359_r8, 0.47359_r8, 0.47359_r8, 0.47359_r8, &
      0.47359_r8, 0.47359_r8, 0.47359_r8, 0.47359_r8, 0.47359_r8, 0.47359_r8, 0.47359_r8, &
      0.47360_r8, 0.47360_r8, 0.47360_r8, 0.47360_r8, 0.47360_r8, 0.47360_r8, 0.47360_r8, &
      0.47360_r8, 0.47360_r8, 0.47360_r8, 0.47360_r8, 0.47360_r8, 0.47360_r8, 0.47361_r8, &
      0.47361_r8, 0.47361_r8, 0.47361_r8, 0.47361_r8, 0.47361_r8, 0.47361_r8, 0.47361_r8, &
      0.47361_r8, 0.47361_r8, 0.47361_r8, 0.47361_r8, 0.47361_r8, 0.47362_r8, 0.47362_r8, &
      0.47362_r8, 0.47362_r8, 0.47362_r8, 0.47362_r8, 0.47362_r8, 0.47362_r8, 0.47362_r8, &
      0.47362_r8, 0.47362_r8, 0.47362_r8, 0.47362_r8, 0.47362_r8, 0.47363_r8, 0.47363_r8, &
      0.47363_r8, 0.47363_r8, 0.47363_r8, 0.47363_r8, 0.47363_r8, 0.47363_r8, 0.47363_r8, &
      0.47363_r8, 0.47363_r8, 0.47363_r8, 0.47363_r8, 0.47364_r8, 0.47364_r8, 0.47364_r8, &
      0.47364_r8, 0.47364_r8, 0.47364_r8, 0.47364_r8, 0.47364_r8, 0.47364_r8, 0.47364_r8, &
      0.47364_r8, 0.47364_r8, 0.47364_r8/)
    real(r8), dimension(11), parameter :: vlimbSidebandFraction27L = &
    (/0.48916_r8, 0.48929_r8, 0.48940_r8, 0.48948_r8, 0.48955_r8, 0.48959_r8, 0.48963_r8, &
      0.48971_r8, 0.48981_r8, 0.48998_r8, 0.49023_r8/)
    real(r8), dimension(4), parameter :: vlimbSidebandFraction32L = &
    (/0.97293_r8, 0.97144_r8, 0.97015_r8, 0.96945_r8/)
    real(r8), dimension(4), parameter :: vlimbSidebandFraction33L = &
    (/0.48936_r8, 0.48045_r8, 0.46610_r8, 0.46746_r8/)
    real(r8), dimension(4), parameter :: vlimbSidebandFraction34L = &
    (/0.95080_r8, 0.94541_r8, 0.94532_r8, 0.94825_r8/)

    real(r8), dimension(25), parameter :: vlimbSidebandFraction2U = &
    (/0.43731_r8, 0.44019_r8, 0.44373_r8, 0.44702_r8, 0.44975_r8, 0.45254_r8, 0.45496_r8, &
      0.45668_r8, 0.45787_r8, 0.45871_r8, 0.45929_r8, 0.45970_r8, 0.45999_r8, 0.46028_r8, &
      0.46069_r8, 0.46126_r8, 0.46205_r8, 0.46316_r8, 0.46469_r8, 0.46675_r8, 0.46896_r8, &
      0.47102_r8, 0.47335_r8, 0.47581_r8, 0.47789_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction3U = &
    (/0.45660_r8, 0.46062_r8, 0.46440_r8, 0.46732_r8, 0.46950_r8, 0.47151_r8, 0.47314_r8, &
      0.47422_r8, 0.47494_r8, 0.47544_r8, 0.47577_r8, 0.47601_r8, 0.47616_r8, 0.47633_r8, &
      0.47656_r8, 0.47686_r8, 0.47729_r8, 0.47786_r8, 0.47861_r8, 0.47958_r8, 0.48052_r8, &
      0.48130_r8, 0.48208_r8, 0.48274_r8, 0.48312_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction4U = &
    (/0.47677_r8, 0.47733_r8, 0.47796_r8, 0.47853_r8, 0.47900_r8, 0.47948_r8, 0.47991_r8, &
      0.48020_r8, 0.48042_r8, 0.48056_r8, 0.48067_r8, 0.48074_r8, 0.48079_r8, 0.48085_r8, &
      0.48092_r8, 0.48102_r8, 0.48116_r8, 0.48136_r8, 0.48162_r8, 0.48198_r8, 0.48235_r8, &
      0.48266_r8, 0.48298_r8, 0.48320_r8, 0.48324_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction5U = &
    (/0.48466_r8, 0.48401_r8, 0.48333_r8, 0.48276_r8, 0.48229_r8, 0.48181_r8, 0.48140_r8, &
      0.48110_r8, 0.48089_r8, 0.48075_r8, 0.48064_r8, 0.48057_r8, 0.48051_r8, 0.48046_r8, &
      0.48039_r8, 0.48029_r8, 0.48014_r8, 0.47993_r8, 0.47964_r8, 0.47924_r8, 0.47880_r8, &
      0.47836_r8, 0.47785_r8, 0.47727_r8, 0.47675_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction6U = &
    (/0.48836_r8, 0.48881_r8, 0.48923_r8, 0.48956_r8, 0.48979_r8, 0.49001_r8, 0.49016_r8, &
      0.49027_r8, 0.49033_r8, 0.49037_r8, 0.49039_r8, 0.49041_r8, 0.49042_r8, 0.49043_r8, &
      0.49045_r8, 0.49047_r8, 0.49049_r8, 0.49051_r8, 0.49052_r8, 0.49049_r8, 0.49038_r8, &
      0.49018_r8, 0.48979_r8, 0.48903_r8, 0.48789_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction7U = &
    (/0.48639_r8, 0.48527_r8, 0.48417_r8, 0.48328_r8, 0.48258_r8, 0.48190_r8, 0.48132_r8, &
      0.48092_r8, 0.48065_r8, 0.48046_r8, 0.48033_r8, 0.48024_r8, 0.48017_r8, 0.48012_r8, &
      0.48002_r8, 0.47990_r8, 0.47973_r8, 0.47949_r8, 0.47918_r8, 0.47881_r8, 0.47844_r8, &
      0.47818_r8, 0.47806_r8, 0.47826_r8, 0.47897_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction8U = &
    (/0.49391_r8, 0.49506_r8, 0.49619_r8, 0.49711_r8, &
      0.49783_r8, 0.49853_r8, 0.49912_r8, 0.49954_r8, 0.49982_r8, &
      0.50001_r8, 0.50016_r8, 0.50025_r8, 0.50032_r8, 0.50039_r8, &
      0.50048_r8, 0.50061_r8, 0.50081_r8, 0.50106_r8, 0.50142_r8, &
      0.50189_r8, 0.50241_r8, 0.50288_r8, 0.50343_r8, 0.50399_r8, &
      0.50446_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction9U = &
    (/0.49235_r8, 0.49301_r8, 0.49375_r8, 0.49437_r8, &
      0.49485_r8, 0.49531_r8, 0.49569_r8, 0.49594_r8, 0.49612_r8, &
      0.49623_r8, 0.49631_r8, 0.49637_r8, 0.49641_r8, 0.49645_r8, &
      0.49651_r8, 0.49658_r8, 0.49669_r8, 0.49684_r8, 0.49703_r8, &
      0.49728_r8, 0.49754_r8, 0.49779_r8, 0.49807_r8, 0.49837_r8, &
      0.49866_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction15U = &
    (/0.04250_r8, 0.04200_r8, 0.04150_r8, 0.04110_r8, &
      0.04080_r8, 0.04040_r8, 0.04010_r8, 0.03990_r8, 0.03980_r8, &
      0.03970_r8, 0.03960_r8, 0.03960_r8, 0.03950_r8, 0.03950_r8, &
      0.03940_r8, 0.03940_r8, 0.03930_r8, 0.03910_r8, 0.03890_r8, &
      0.03860_r8, 0.03830_r8, 0.03790_r8, 0.03750_r8, 0.03700_r8, &
      0.03650_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction16U = &
    (/0.19600_r8, 0.19130_r8, 0.18710_r8, 0.18330_r8, &
      0.18040_r8, 0.17730_r8, 0.17480_r8, 0.17290_r8, 0.17160_r8, &
      0.17060_r8, 0.17000_r8, 0.16950_r8, 0.16920_r8, 0.16890_r8, &
      0.16840_r8, 0.16780_r8, 0.16690_r8, 0.16560_r8, 0.16370_r8, &
      0.16110_r8, 0.15830_r8, 0.15510_r8, 0.15160_r8, 0.14680_r8, &
      0.14280_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction17U = &
    (/0.57500_r8, 0.58080_r8, 0.58650_r8, 0.59150_r8, &
      0.59550_r8, 0.59970_r8, 0.60320_r8, 0.60560_r8, 0.60730_r8, &
      0.60850_r8, 0.60940_r8, 0.61000_r8, 0.61050_r8, 0.61090_r8, &
      0.61150_r8, 0.61240_r8, 0.61360_r8, 0.61540_r8, 0.61780_r8, &
      0.62140_r8, 0.62510_r8, 0.62910_r8, 0.63410_r8, 0.63990_r8, &
      0.64600_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction18U = &
    (/0.05350_r8, 0.05320_r8, 0.05300_r8, 0.05270_r8, &
      0.05260_r8, 0.05240_r8, 0.05220_r8, 0.05210_r8, 0.05200_r8, &
      0.05200_r8, 0.05200_r8, 0.05190_r8, 0.05190_r8, 0.05190_r8, &
      0.05190_r8, 0.05180_r8, 0.05180_r8, 0.05170_r8, 0.05160_r8, &
      0.05140_r8, 0.05130_r8, 0.05110_r8, 0.05090_r8, 0.05060_r8, &
      0.05030_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction19U = &
    (/0.21500_r8, 0.21020_r8, 0.20590_r8, 0.20200_r8, &
      0.19890_r8, 0.19600_r8, 0.19330_r8, 0.19140_r8, 0.19010_r8, &
      0.18910_r8, 0.18840_r8, 0.18800_r8, 0.18760_r8, 0.18730_r8, &
      0.18680_r8, 0.18620_r8, 0.18530_r8, 0.18400_r8, 0.18190_r8, &
      0.17940_r8, 0.17630_r8, 0.17330_r8, 0.16960_r8, 0.16490_r8, &
      0.16040_r8/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction20U = &
    (/0.58130_r8, 0.58690_r8, 0.59350_r8, 0.59800_r8, &
      0.60200_r8, 0.60630_r8, 0.60960_r8, 0.61210_r8, 0.61380_r8, &
      0.61500_r8, 0.61590_r8, 0.61650_r8, 0.61700_r8, 0.61740_r8, &
      0.61800_r8, 0.61880_r8, 0.62010_r8, 0.62190_r8, 0.62430_r8, &
      0.62800_r8, 0.63150_r8, 0.63550_r8, 0.64070_r8, 0.64620_r8, &
      0.65240_r8/)
    real(r8), dimension(129), parameter :: vlimbSidebandFraction23U = &
    (/0.46025_r8, 0.46024_r8, 0.46024_r8, 0.46024_r8, &
      0.46023_r8, 0.46023_r8, 0.46023_r8, 0.46022_r8, 0.46022_r8, &
      0.46022_r8, 0.46021_r8, 0.46021_r8, 0.46021_r8, 0.46020_r8, &
      0.46020_r8, 0.46020_r8, 0.46019_r8, 0.46019_r8, 0.46019_r8, &
      0.46018_r8, 0.46018_r8, 0.46018_r8, 0.46017_r8, 0.46017_r8, &
      0.46017_r8, 0.46016_r8, 0.46016_r8, 0.46016_r8, 0.46015_r8, &
      0.46015_r8, 0.46015_r8, 0.46014_r8, 0.46014_r8, 0.46014_r8, &
      0.46013_r8, 0.46013_r8, 0.46013_r8, 0.46012_r8, 0.46012_r8, &
      0.46012_r8, 0.46011_r8, 0.46011_r8, 0.46011_r8, 0.46010_r8, &
      0.46010_r8, 0.46010_r8, 0.46009_r8, 0.46009_r8, 0.46009_r8, &
      0.46008_r8, 0.46008_r8, 0.46008_r8, 0.46007_r8, 0.46007_r8, &
      0.46007_r8, 0.46006_r8, 0.46006_r8, 0.46006_r8, 0.46005_r8, &
      0.46005_r8, 0.46005_r8, 0.46004_r8, 0.46004_r8, 0.46004_r8, &
      0.46003_r8, 0.46003_r8, 0.46003_r8, 0.46002_r8, 0.46002_r8, &
      0.46002_r8, 0.46001_r8, 0.46001_r8, 0.46001_r8, 0.46000_r8, &
      0.46000_r8, 0.46000_r8, 0.45999_r8, 0.45999_r8, 0.45999_r8, &
      0.45998_r8, 0.45998_r8, 0.45998_r8, 0.45997_r8, 0.45997_r8, &
      0.45997_r8, 0.45996_r8, 0.45996_r8, 0.45996_r8, 0.45995_r8, &
      0.45995_r8, 0.45995_r8, 0.45994_r8, 0.45994_r8, 0.45994_r8, &
      0.45994_r8, 0.45994_r8, 0.45994_r8, 0.45993_r8, 0.45993_r8, &
      0.45993_r8, 0.45992_r8, 0.45992_r8, 0.45991_r8, 0.45991_r8, &
      0.45991_r8, 0.45990_r8, 0.45990_r8, 0.45990_r8, 0.45989_r8, &
      0.45989_r8, 0.45989_r8, 0.45988_r8, 0.45988_r8, 0.45988_r8, &
      0.45987_r8, 0.45987_r8, 0.45987_r8, 0.45986_r8, 0.45986_r8, &
      0.45986_r8, 0.45985_r8, 0.45985_r8, 0.45985_r8, 0.45984_r8, &
      0.45984_r8, 0.45984_r8, 0.45983_r8, 0.45983_r8, 0.45983_r8/)
    real(r8), dimension(129), parameter :: vlimbSidebandFraction24U = &
    (/0.48014_r8, 0.48014_r8, 0.48014_r8, 0.48014_r8, &
      0.48014_r8, 0.48014_r8, 0.48014_r8, 0.48014_r8, 0.48014_r8, &
      0.48014_r8, 0.48014_r8, 0.48014_r8, 0.48014_r8, 0.48014_r8, &
      0.48014_r8, 0.48014_r8, 0.48014_r8, 0.48015_r8, 0.48015_r8, &
      0.48015_r8, 0.48015_r8, 0.48015_r8, 0.48015_r8, 0.48015_r8, &
      0.48015_r8, 0.48015_r8, 0.48015_r8, 0.48015_r8, 0.48015_r8, &
      0.48015_r8, 0.48015_r8, 0.48015_r8, 0.48015_r8, 0.48015_r8, &
      0.48015_r8, 0.48015_r8, 0.48015_r8, 0.48015_r8, 0.48015_r8, &
      0.48015_r8, 0.48016_r8, 0.48016_r8, 0.48016_r8, 0.48016_r8, &
      0.48016_r8, 0.48016_r8, 0.48016_r8, 0.48016_r8, 0.48016_r8, &
      0.48016_r8, 0.48016_r8, 0.48016_r8, 0.48016_r8, 0.48016_r8, &
      0.48016_r8, 0.48016_r8, 0.48016_r8, 0.48016_r8, 0.48016_r8, &
      0.48016_r8, 0.48016_r8, 0.48016_r8, 0.48016_r8, 0.48017_r8, &
      0.48017_r8, 0.48017_r8, 0.48017_r8, 0.48017_r8, 0.48017_r8, &
      0.48017_r8, 0.48017_r8, 0.48017_r8, 0.48017_r8, 0.48017_r8, &
      0.48017_r8, 0.48017_r8, 0.48017_r8, 0.48017_r8, 0.48017_r8, &
      0.48017_r8, 0.48017_r8, 0.48017_r8, 0.48017_r8, 0.48017_r8, &
      0.48017_r8, 0.48017_r8, 0.48018_r8, 0.48018_r8, 0.48018_r8, &
      0.48018_r8, 0.48018_r8, 0.48018_r8, 0.48018_r8, 0.48018_r8, &
      0.48018_r8, 0.48018_r8, 0.48018_r8, 0.48018_r8, 0.48018_r8, &
      0.48018_r8, 0.48018_r8, 0.48018_r8, 0.48018_r8, 0.48018_r8, &
      0.48018_r8, 0.48018_r8, 0.48018_r8, 0.48018_r8, 0.48018_r8, &
      0.48019_r8, 0.48019_r8, 0.48019_r8, 0.48019_r8, 0.48019_r8, &
      0.48019_r8, 0.48019_r8, 0.48019_r8, 0.48019_r8, 0.48019_r8, &
      0.48019_r8, 0.48019_r8, 0.48019_r8, 0.48019_r8, 0.48019_r8, &
      0.48019_r8, 0.48019_r8, 0.48019_r8, 0.48019_r8, 0.48019_r8/)
    real(r8), dimension(129), parameter :: vlimbSidebandFraction25U = &
    (/0.49647_r8, 0.49647_r8, 0.49647_r8, 0.49647_r8, &
      0.49647_r8, 0.49647_r8, 0.49647_r8, 0.49647_r8, 0.49647_r8, &
      0.49647_r8, 0.49646_r8, 0.49646_r8, 0.49646_r8, 0.49646_r8, &
      0.49646_r8, 0.49646_r8, 0.49646_r8, 0.49646_r8, 0.49646_r8, &
      0.49646_r8, 0.49646_r8, 0.49646_r8, 0.49646_r8, 0.49646_r8, &
      0.49645_r8, 0.49645_r8, 0.49645_r8, 0.49645_r8, 0.49645_r8, &
      0.49645_r8, 0.49645_r8, 0.49645_r8, 0.49645_r8, 0.49645_r8, &
      0.49645_r8, 0.49645_r8, 0.49645_r8, 0.49644_r8, 0.49644_r8, &
      0.49644_r8, 0.49644_r8, 0.49644_r8, 0.49644_r8, 0.49644_r8, &
      0.49644_r8, 0.49644_r8, 0.49644_r8, 0.49644_r8, 0.49644_r8, &
      0.49644_r8, 0.49643_r8, 0.49643_r8, 0.49643_r8, 0.49643_r8, &
      0.49643_r8, 0.49643_r8, 0.49643_r8, 0.49643_r8, 0.49643_r8, &
      0.49643_r8, 0.49643_r8, 0.49643_r8, 0.49643_r8, 0.49642_r8, &
      0.49642_r8, 0.49642_r8, 0.49642_r8, 0.49642_r8, 0.49642_r8, &
      0.49642_r8, 0.49642_r8, 0.49642_r8, 0.49642_r8, 0.49642_r8, &
      0.49642_r8, 0.49642_r8, 0.49641_r8, 0.49641_r8, 0.49641_r8, &
      0.49641_r8, 0.49641_r8, 0.49641_r8, 0.49641_r8, 0.49641_r8, &
      0.49641_r8, 0.49641_r8, 0.49641_r8, 0.49641_r8, 0.49641_r8, &
      0.49640_r8, 0.49640_r8, 0.49640_r8, 0.49640_r8, 0.49640_r8, &
      0.49640_r8, 0.49640_r8, 0.49640_r8, 0.49640_r8, 0.49640_r8, &
      0.49640_r8, 0.49640_r8, 0.49640_r8, 0.49640_r8, 0.49639_r8, &
      0.49639_r8, 0.49639_r8, 0.49639_r8, 0.49639_r8, 0.49639_r8, &
      0.49639_r8, 0.49639_r8, 0.49639_r8, 0.49639_r8, 0.49639_r8, &
      0.49639_r8, 0.49639_r8, 0.49638_r8, 0.49638_r8, 0.49638_r8, &
      0.49638_r8, 0.49638_r8, 0.49638_r8, 0.49638_r8, 0.49638_r8, &
      0.49638_r8, 0.49638_r8, 0.49638_r8, 0.49638_r8, 0.49638_r8/)
    real(r8), dimension(11), parameter :: vlimbSidebandFraction27U = &
    (/0.49004_r8, 0.48991_r8, 0.48981_r8, 0.48973_r8, &
      0.48966_r8, 0.48962_r8, 0.48957_r8, 0.48949_r8, 0.48939_r8, &
      0.48923_r8, 0.48897_r8/)
    real(r8), dimension(4), parameter :: vlimbSidebandFraction33U = &
    (/0.48061_r8, 0.48947_r8, 0.50383_r8, 0.50247_r8/)

    integer :: limbSidebandFraction1L = 0
    integer :: limbSidebandFraction2L = 0
    integer :: limbSidebandFraction3L = 0
    integer :: limbSidebandFraction4L = 0
    integer :: limbSidebandFraction5L = 0
    integer :: limbSidebandFraction6L = 0
    integer :: limbSidebandFraction7L = 0
    integer :: limbSidebandFraction8L = 0
    integer :: limbSidebandFraction9L = 0
    integer :: limbSidebandFraction10L = 0
    integer :: limbSidebandFraction11L = 0
    integer :: limbSidebandFraction12L = 0
    integer :: limbSidebandFraction13L = 0
    integer :: limbSidebandFraction14L = 0
    integer :: limbSidebandFraction15L = 0
    integer :: limbSidebandFraction16L = 0
    integer :: limbSidebandFraction17L = 0
    integer :: limbSidebandFraction18L = 0
    integer :: limbSidebandFraction19L = 0
    integer :: limbSidebandFraction20L = 0
    integer :: limbSidebandFraction21L = 0
    integer :: limbSidebandFraction22L = 0
    integer :: limbSidebandFraction23L = 0
    integer :: limbSidebandFraction24L = 0
    integer :: limbSidebandFraction25L = 0
    integer :: limbSidebandFraction26L = 0
    integer :: limbSidebandFraction27L = 0
    integer :: limbSidebandFraction28L = 0
    integer :: limbSidebandFraction29L = 0
    integer :: limbSidebandFraction30L = 0
    integer :: limbSidebandFraction31L = 0
    integer :: limbSidebandFraction32L = 0
    integer :: limbSidebandFraction33L = 0
    integer :: limbSidebandFraction34L = 0

    integer :: limbSidebandFraction2U = 0
    integer :: limbSidebandFraction3U = 0
    integer :: limbSidebandFraction4U = 0
    integer :: limbSidebandFraction5U = 0
    integer :: limbSidebandFraction6U = 0
    integer :: limbSidebandFraction7U = 0
    integer :: limbSidebandFraction8U = 0
    integer :: limbSidebandFraction9U = 0
    integer :: limbSidebandFraction10U = 0
    integer :: limbSidebandFraction11U = 0
    integer :: limbSidebandFraction12U = 0
    integer :: limbSidebandFraction13U = 0
    integer :: limbSidebandFraction14U = 0
    integer :: limbSidebandFraction15U = 0
    integer :: limbSidebandFraction16U = 0
    integer :: limbSidebandFraction17U = 0
    integer :: limbSidebandFraction18U = 0
    integer :: limbSidebandFraction19U = 0
    integer :: limbSidebandFraction20U = 0
    integer :: limbSidebandFraction23U = 0
    integer :: limbSidebandFraction24U = 0
    integer :: limbSidebandFraction25U = 0
    integer :: limbSidebandFraction27U = 0
    integer :: limbSidebandFraction28U = 0
    integer :: limbSidebandFraction29U = 0
    integer :: limbSidebandFraction30U = 0
    integer :: limbSidebandFraction31U = 0
    integer :: limbSidebandFraction33U = 0

    contains

    ! InitQuantityTemplates should be called priori to this subroutine
    subroutine CreateLimbSidebandFractions (fakeChunk, filedatabase, qtyTemplates)
       use CFM_QuantityTemplate_m, only: CreateQtyTemplate, AddQuantityTemplateToDatabase, &
                                         QuantityTemplate_T
       use INIT_TABLES_MODULE, only: l_limbsidebandFraction
       use Chunks_m, only: MLSChunk_T
       use MLSCommon, only: MLSFile_T

       type(MLSChunk_T), intent(in) :: fakeChunk
       type (MLSFile_T), dimension(:), pointer :: filedatabase
       type (QuantityTemplate_T), dimension(:), pointer :: qtyTemplates

       type(QuantityTemplate_T) :: lsbFraction

       ! Executables

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R1A:118.B1LF:PT")
       limbSidebandFraction1L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R2:190.B2LF:H2O")
       limbSidebandFraction2L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R2:190.B3LF:N2O")
       limbSidebandFraction3L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R2:190.B4LF:HNO3")
       limbSidebandFraction4L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R2:190.B5LF:ClO")
       limbSidebandFraction5L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R2:190.B6LF:O3")
       limbSidebandFraction6L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R3:240.B7LF:O3")
       limbSidebandFraction7L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R3:240.B8LF:PT")
       limbSidebandFraction8L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R3:240.B9LF:CO")
       limbSidebandFraction9L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R4:640.B10LF:ClO")
       limbSidebandFraction10L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R4:640.B11LF:BrO")
       limbSidebandFraction11L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R4:640.B12LF:N2O")
       limbSidebandFraction12L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R4:640.B13LF:HCl")
       limbSidebandFraction13L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R4:640.B14LF:O3")
       limbSidebandFraction14L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R5H:2T5.B15LF:OH")
       limbSidebandFraction15L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R5H:2T5.B16LF:OH")
       limbSidebandFraction16L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R5H:2T5.B17LF:PT")
       limbSidebandFraction17L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R5V:2T5.B18LF:OH")
       limbSidebandFraction18L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R5V:2T5.B19LF:OH")
       limbSidebandFraction19L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R5V:2T5.B20LF:PT")
       limbSidebandFraction20L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R1B:118.B21LF:PT")
       limbSidebandFraction21L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R1A:118.B22LD:PT")
       limbSidebandFraction22L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R2:190.B23LD:H2O")
       limbSidebandFraction23L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R3:240.B24LD:O3")
       limbSidebandFraction24L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R3:240.B25LD:CO")
       limbSidebandFraction25L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R1B:118.B26LD:PT")
       limbSidebandFraction26L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R2:190.B27LM:HCN")
       limbSidebandFraction27L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R4:640.B28LM:HO2")
       limbSidebandFraction28L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R4:640.B29LM:HOCl")
       limbSidebandFraction29L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R4:640.B30LM:HO2")
       limbSidebandFraction30L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R4:640.B31LM:BrO")
       limbSidebandFraction31L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R1A:118.B32LW:PT")
       limbSidebandFraction32L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R3:240.B33LW:O3")
       limbSidebandFraction33L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R1B:118.B34LW:PT")
       limbSidebandFraction34L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       ! Now the the upper band
       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R2:190.B2UF:H2O")
       limbSidebandFraction2U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R2:190.B3UF:N2O")
       limbSidebandFraction3U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R2:190.B4UF:HNO3")
       limbSidebandFraction4U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R2:190.B5UF:ClO")
       limbSidebandFraction5U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R2:190.B6UF:O3")
       limbSidebandFraction6U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R3:240.B7UF:O3")
       limbSidebandFraction7U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R3:240.B8UF:PT")
       limbSidebandFraction8U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R3:240.B9UF:CO")
       limbSidebandFraction9U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R4:640.B10UF:ClO")
       limbSidebandFraction10U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R4:640.B11UF:BrO")
       limbSidebandFraction11U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R4:640.B12UF:N2O")
       limbSidebandFraction12U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R4:640.B13UF:HCl")
       limbSidebandFraction13U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R4:640.B14UF:O3")
       limbSidebandFraction14U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R5H:2T5.B15UF:OH")
       limbSidebandFraction15U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R5H:2T5.B16UF:OH")
       limbSidebandFraction16U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R5H:2T5.B17UF:PT")
       limbSidebandFraction17U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R5V:2T5.B18UF:OH")
       limbSidebandFraction18U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R5V:2T5.B19UF:OH")
       limbSidebandFraction19U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R5V:2T5.B20UF:PT")
       limbSidebandFraction20U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R2:190.B23UD:H2O")
       limbSidebandFraction23U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R3:240.B24UD:O3")
       limbSidebandFraction24U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R3:240.B25UD:CO")
       limbSidebandFraction25U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R2:190.B27UM:HCN")
       limbSidebandFraction27U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R4:640.B28UM:HO2")
       limbSidebandFraction28U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R4:640.B29UM:HOCl")
       limbSidebandFraction29U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R4:640.B30UM:HO2")
       limbSidebandFraction30U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R4:640.B31UM:BrO")
       limbSidebandFraction31U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=fakeChunk, qSignal="R3:240.B33UW:O3")
       limbSidebandFraction33U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

    end subroutine

    subroutine FillLimbSidebandFractions (stateVectorExtra)
       use CFM_Vector_m, only: Vector_T, VectorValue_T, Dump
       use CFM_Vector_m, only: GetVectorQtyByTemplateIndex
       use CFM_Fill_M, only: ExplicitFillVectorQuantity, SpreadFillVectorQuantity

       type (Vector_T), intent(in) :: stateVectorExtra
       type(VectorValue_T) :: quantity

       ! Executables
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction1L)
       call SpreadFillVectorQuantity (quantity, 0.97333_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction2L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction2L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction3L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction3L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction7L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction7L)
       !print *, "LimbSidebandFraction7L value"
       !call dump(quantity, details=3)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction8L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction8L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction9L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction9L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction10L)
       call SpreadFillVectorQuantity (quantity, 0.50882_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction11L)
       call SpreadFillVectorQuantity (quantity, 0.50882_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction12L)
       call SpreadFillVectorQuantity (quantity, 0.50690_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction13L)
       call SpreadFillVectorQuantity (quantity, 0.51563_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction14L)
       call SpreadFillVectorQuantity (quantity, 0.51563_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction15L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction15L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction16L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction16L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction17L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction17L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction18L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction18L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction19L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction19L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction20L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction20L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction21L)
       call SpreadFillVectorQuantity (quantity, 0.94708_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction22L)
       call SpreadFillVectorQuantity (quantity, 0.97333_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction23L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction23L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction24L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction24L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction25L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction25L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction26L)
       call SpreadFillVectorQuantity (quantity, 0.94708_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction27L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction27L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction28L)
       call SpreadFillVectorQuantity (quantity, 0.50882_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction29L)
       call SpreadFillVectorQuantity (quantity, 0.50882_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction30L)
       call SpreadFillVectorQuantity (quantity, 0.51563_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction31L)
       call SpreadFillVectorQuantity (quantity, 0.51563_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction32L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction32L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction33L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction33L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction34L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction34L)

       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction2U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction2U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction3U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction3U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction4U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction4U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction5U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction5U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction6U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction6U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction7U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction7U)
       !print *, "LimbSidebandFraction7U value"
       !call dump(quantity, details=3)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction8U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction8U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction9U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction9U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction10U)
       call SpreadFillVectorQuantity (quantity, 0.46045_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction11U)
       call SpreadFillVectorQuantity (quantity, 0.46045_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction12U)
       call SpreadFillVectorQuantity (quantity, 0.46232_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction13U)
       call SpreadFillVectorQuantity (quantity, 0.45360_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction14U)
       call SpreadFillVectorQuantity (quantity, 0.45360_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction15U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction15U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction16U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction16U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction17U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction17U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction18U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction18U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction19U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction19U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction20U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction20U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction23U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction23U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction24U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction24U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction25U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction25U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction27U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction27U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction28U)
       call SpreadFillVectorQuantity (quantity, 0.46045_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction29U)
       call SpreadFillVectorQuantity (quantity, 0.46045_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction30U)
       call SpreadFillVectorQuantity (quantity, 0.45360_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction31U)
       call SpreadFillVectorQuantity (quantity, 0.45360_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction33U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction33U)

    end subroutine

!--------------------------- end bloc --------------------------------------
   logical function not_used_here()
   character (len=*), parameter :: IdParm = &
       "$Id$"
   character (len=len(idParm)) :: Id = idParm
      not_used_here = (id(1:1) == ModuleName(1:1))
      print *, Id ! .mod files sometimes change if PRINT is added
   end function not_used_here
!---------------------------------------------------------------------------
end module
