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
    (/0.54213, 0.53927, 0.53572, 0.53244, 0.52971, 0.52693, 0.52451, &
      0.52279, 0.52160, 0.52076, 0.52018, 0.51977, 0.51948, 0.51919, &
      0.51878, 0.51822, 0.51742, 0.51632, 0.51479, 0.51273, 0.51052, &
      0.50846, 0.50613, 0.50368, 0.50160/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction3L = &
    (/0.52287, 0.51885, 0.51507, 0.51216, 0.50998, 0.50797, 0.50635, &
      0.50526, 0.50455, 0.50405, 0.50372, 0.50348, 0.50332, 0.50316, &
      0.50293, 0.50263, 0.50220, 0.50163, 0.50088, 0.49992, 0.49898, &
      0.49819, 0.49741, 0.49676, 0.49637/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction7L = &
    (/0.48358, 0.48470, 0.48580, 0.48669, 0.48739, 0.48807, 0.48865, &
      0.48905, 0.48932, 0.48950, 0.48964, 0.48973, 0.48979, 0.48985, &
      0.48995, 0.49007, 0.49024, 0.49047, 0.49078, 0.49116, 0.49153, &
      0.49178, 0.49191, 0.49171, 0.49100/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction8L = &
    (/0.47597, 0.47481, 0.47369, 0.47277, &
      0.47205, 0.47135, 0.47076, 0.47034, 0.47006, &
      0.46987, 0.46972, 0.46962, 0.46956, 0.46949, &
      0.46940, 0.46927, 0.46907, 0.46882, 0.46846, &
      0.46799, 0.46747, 0.46700, 0.46645, 0.46589, &
      0.46543/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction9L = &
    (/0.47767, 0.47701, 0.47627, 0.47565, &
      0.47517, 0.47471, 0.47433, 0.47408, 0.47390, &
      0.47379, 0.47371, 0.47365, 0.47361, 0.47357, &
      0.47352, 0.47344, 0.47333, 0.47319, 0.47299, &
      0.47274, 0.47248, 0.47224, 0.47195, 0.47165, &
      0.47136/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction15L = &
    (/0.95750, 0.95800, 0.95850, 0.95890, &
      0.95920, 0.95960, 0.95990, 0.96010, 0.96020, &
      0.96030, 0.96040, 0.96040, 0.96050, 0.96050, &
      0.96060, 0.96060, 0.96070, 0.96090, 0.96110, &
      0.96140, 0.96170, 0.96210, 0.96250, 0.96300, &
      0.96350/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction16L = &
    (/0.80400, 0.80870, 0.81290, 0.81670, &
      0.81960, 0.82270, 0.82520, 0.82710, 0.82840, &
      0.82940, 0.83000, 0.83050, 0.83080, 0.83110, &
      0.83160, 0.83220, 0.83310, 0.83440, 0.83630, &
      0.83890, 0.84170, 0.84490, 0.84840, 0.85320, &
      0.85720/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction17L = &
    (/0.42500, 0.41920, 0.41350, 0.40850, &
      0.40450, 0.40030, 0.39680, 0.39440, 0.39270, &
      0.39150, 0.39060, 0.39000, 0.38950, 0.38910, &
      0.38850, 0.38760, 0.38640, 0.38460, 0.38220, &
      0.37860, 0.37490, 0.37090, 0.36590, 0.36010, &
      0.35400/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction18L = &
    (/0.94650, 0.94680, 0.94700, 0.94730, &
      0.94740, 0.94760, 0.94780, 0.94790, 0.94800, &
      0.94800, 0.94800, 0.94810, 0.94810, 0.94810, &
      0.94810, 0.94820, 0.94820, 0.94830, 0.94840, &
      0.94860, 0.94870, 0.94890, 0.94910, 0.94940, &
      0.94970/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction19L = &
    (/0.78500, 0.78980, 0.79410, 0.79800, &
      0.80110, 0.80400, 0.80670, 0.80860, 0.80990, &
      0.81090, 0.81160, 0.81200, 0.81240, 0.81270, &
      0.81320, 0.81380, 0.81470, 0.81600, 0.81810, &
      0.82060, 0.82370, 0.82670, 0.83040, 0.83510, &
      0.83960/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction20L = &
    (/0.41870, 0.41310, 0.40650, 0.40200, &
      0.39800, 0.39370, 0.39040, 0.38790, 0.38620, &
      0.38500, 0.38410, 0.38350, 0.38300, 0.38260, &
      0.38200, 0.38120, 0.37990, 0.37810, 0.37570, &
      0.37200, 0.36850, 0.36450, 0.35930, 0.35380, &
      0.34760/)
    real(r8), dimension(129), parameter :: vlimbSidebandFraction23L = &
    (/0.51922, 0.51923, 0.51923, 0.51923, &
      0.51924, 0.51924, 0.51924, 0.51925, 0.51925, &
      0.51925, 0.51926, 0.51926, 0.51926, 0.51927, &
      0.51927, 0.51927, 0.51928, 0.51928, 0.51928, &
      0.51929, 0.51929, 0.51929, 0.51930, 0.51930, &
      0.51930, 0.51931, 0.51931, 0.51931, 0.51932, &
      0.51932, 0.51932, 0.51933, 0.51933, 0.51933, &
      0.51934, 0.51934, 0.51934, 0.51935, 0.51935, &
      0.51935, 0.51936, 0.51936, 0.51936, 0.51937, &
      0.51937, 0.51937, 0.51938, 0.51938, 0.51938, &
      0.51939, 0.51939, 0.51939, 0.51940, 0.51940, &
      0.51940, 0.51941, 0.51941, 0.51941, 0.51942, &
      0.51942, 0.51942, 0.51943, 0.51943, 0.51943, &
      0.51944, 0.51944, 0.51944, 0.51945, 0.51945, &
      0.51945, 0.51946, 0.51946, 0.51946, 0.51947, &
      0.51947, 0.51947, 0.51948, 0.51948, 0.51948, &
      0.51949, 0.51949, 0.51949, 0.51950, 0.51950, &
      0.51950, 0.51951, 0.51951, 0.51951, 0.51952, &
      0.51952, 0.51952, 0.51953, 0.51953, 0.51953, &
      0.51954, 0.51954, 0.51954, 0.51955, 0.51955, &
      0.51955, 0.51956, 0.51956, 0.51957, 0.51957, &
      0.51957, 0.51958, 0.51958, 0.51958, 0.51959, &
      0.51959, 0.51959, 0.51960, 0.51960, 0.51960, &
      0.51961, 0.51961, 0.51961, 0.51962, 0.51962, &
      0.51962, 0.51963, 0.51963, 0.51963, 0.51964, &
      0.51964, 0.51964, 0.51965, 0.51965, 0.51965/)
    real(r8), dimension(129), parameter :: vlimbSidebandFraction24L = &
    (/0.48982, 0.48982, 0.48982, 0.48982, &
      0.48982, 0.48982, 0.48982, 0.48982, 0.48982, &
      0.48982, 0.48982, 0.48982, 0.48982, 0.48982, &
      0.48982, 0.48982, 0.48982, 0.48981, 0.48981, &
      0.48981, 0.48981, 0.48981, 0.48981, 0.48981, &
      0.48981, 0.48981, 0.48981, 0.48981, 0.48981, &
      0.48981, 0.48981, 0.48981, 0.48981, 0.48981, &
      0.48981, 0.48981, 0.48981, 0.48981, 0.48981, &
      0.48981, 0.48980, 0.48980, 0.48980, 0.48980, &
      0.48980, 0.48980, 0.48980, 0.48980, 0.48980, &
      0.48980, 0.48980, 0.48980, 0.48980, 0.48980, &
      0.48980, 0.48980, 0.48980, 0.48980, 0.48980, &
      0.48980, 0.48980, 0.48980, 0.48980, 0.48979, &
      0.48979, 0.48979, 0.48979, 0.48979, 0.48979, &
      0.48979, 0.48979, 0.48979, 0.48979, 0.48979, &
      0.48979, 0.48979, 0.48979, 0.48979, 0.48979, &
      0.48979, 0.48979, 0.48979, 0.48979, 0.48979, &
      0.48979, 0.48979, 0.48979, 0.48979, 0.48979, &
      0.48979, 0.48979, 0.48979, 0.48979, 0.48979, &
      0.48979, 0.48979, 0.48979, 0.48979, 0.48979, &
      0.48979, 0.48979, 0.48979, 0.48979, 0.48979, &
      0.48979, 0.48979, 0.48979, 0.48979, 0.48979, &
      0.48978, 0.48978, 0.48978, 0.48978, 0.48978, &
      0.48978, 0.48978, 0.48978, 0.48978, 0.48978, &
      0.48978, 0.48978, 0.48978, 0.48978, 0.48978, &
      0.48978, 0.48978, 0.48978, 0.48978, 0.48978/)
    real(r8), dimension(129), parameter :: vlimbSidebandFraction25L = &
    (/0.47355, 0.47355, 0.47355, 0.47355, 0.47355, 0.47355, 0.47355, &
      0.47355, 0.47355, 0.47355, 0.47356, 0.47356, 0.47356, 0.47356, &
      0.47356, 0.47356, 0.47356, 0.47356, 0.47356, 0.47356, 0.47356, &
      0.47356, 0.47356, 0.47356, 0.47357, 0.47357, 0.47357, 0.47357, &
      0.47357, 0.47357, 0.47357, 0.47357, 0.47357, 0.47357, 0.47357, &
      0.47357, 0.47357, 0.47358, 0.47358, 0.47358, 0.47358, 0.47358, &
      0.47358, 0.47358, 0.47358, 0.47358, 0.47358, 0.47358, 0.47358, &
      0.47358, 0.47359, 0.47359, 0.47359, 0.47359, 0.47359, 0.47359, &
      0.47359, 0.47359, 0.47359, 0.47359, 0.47359, 0.47359, 0.47359, &
      0.47360, 0.47360, 0.47360, 0.47360, 0.47360, 0.47360, 0.47360, &
      0.47360, 0.47360, 0.47360, 0.47360, 0.47360, 0.47360, 0.47361, &
      0.47361, 0.47361, 0.47361, 0.47361, 0.47361, 0.47361, 0.47361, &
      0.47361, 0.47361, 0.47361, 0.47361, 0.47361, 0.47362, 0.47362, &
      0.47362, 0.47362, 0.47362, 0.47362, 0.47362, 0.47362, 0.47362, &
      0.47362, 0.47362, 0.47362, 0.47362, 0.47362, 0.47363, 0.47363, &
      0.47363, 0.47363, 0.47363, 0.47363, 0.47363, 0.47363, 0.47363, &
      0.47363, 0.47363, 0.47363, 0.47363, 0.47364, 0.47364, 0.47364, &
      0.47364, 0.47364, 0.47364, 0.47364, 0.47364, 0.47364, 0.47364, &
      0.47364, 0.47364, 0.47364/)
    real(r8), dimension(11), parameter :: vlimbSidebandFraction27L = &
    (/0.48916, 0.48929, 0.48940, 0.48948, 0.48955, 0.48959, 0.48963, &
      0.48971, 0.48981, 0.48998, 0.49023/)
    real(r8), dimension(4), parameter :: vlimbSidebandFraction32L = &
    (/0.97293, 0.97144, 0.97015, 0.96945/)
    real(r8), dimension(4), parameter :: vlimbSidebandFraction33L = &
    (/0.48936, 0.48045, 0.46610, 0.46746/)
    real(r8), dimension(4), parameter :: vlimbSidebandFraction34L = &
    (/0.95080, 0.94541, 0.94532, 0.94825/)

    real(r8), dimension(25), parameter :: vlimbSidebandFraction2U = &
    (/0.43731, 0.44019, 0.44373, 0.44702, 0.44975, 0.45254, 0.45496, &
      0.45668, 0.45787, 0.45871, 0.45929, 0.45970, 0.45999, 0.46028, &
      0.46069, 0.46126, 0.46205, 0.46316, 0.46469, 0.46675, 0.46896, &
      0.47102, 0.47335, 0.47581, 0.47789/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction3U = &
    (/0.45660, 0.46062, 0.46440, 0.46732, 0.46950, 0.47151, 0.47314, &
      0.47422, 0.47494, 0.47544, 0.47577, 0.47601, 0.47616, 0.47633, &
      0.47656, 0.47686, 0.47729, 0.47786, 0.47861, 0.47958, 0.48052, &
      0.48130, 0.48208, 0.48274, 0.48312/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction4U = &
    (/0.47677, 0.47733, 0.47796, 0.47853, 0.47900, 0.47948, 0.47991, &
      0.48020, 0.48042, 0.48056, 0.48067, 0.48074, 0.48079, 0.48085, &
      0.48092, 0.48102, 0.48116, 0.48136, 0.48162, 0.48198, 0.48235, &
      0.48266, 0.48298, 0.48320, 0.48324/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction5U = &
    (/0.48466, 0.48401, 0.48333, 0.48276, 0.48229, 0.48181, 0.48140, &
      0.48110, 0.48089, 0.48075, 0.48064, 0.48057, 0.48051, 0.48046, &
      0.48039, 0.48029, 0.48014, 0.47993, 0.47964, 0.47924, 0.47880, &
      0.47836, 0.47785, 0.47727, 0.47675/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction6U = &
    (/0.48836, 0.48881, 0.48923, 0.48956, 0.48979, 0.49001, 0.49016, &
      0.49027, 0.49033, 0.49037, 0.49039, 0.49041, 0.49042, 0.49043, &
      0.49045, 0.49047, 0.49049, 0.49051, 0.49052, 0.49049, 0.49038, &
      0.49018, 0.48979, 0.48903, 0.48789/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction7U = &
    (/0.48639, 0.48527, 0.48417, 0.48328, 0.48258, 0.48190, 0.48132, &
      0.48092, 0.48065, 0.48046, 0.48033, 0.48024, 0.48017, 0.48012, &
      0.48002, 0.47990, 0.47973, 0.47949, 0.47918, 0.47881, 0.47844, &
      0.47818, 0.47806, 0.47826, 0.47897/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction8U = &
    (/0.49391, 0.49506, 0.49619, 0.49711, &
      0.49783, 0.49853, 0.49912, 0.49954, 0.49982, &
      0.50001, 0.50016, 0.50025, 0.50032, 0.50039, &
      0.50048, 0.50061, 0.50081, 0.50106, 0.50142, &
      0.50189, 0.50241, 0.50288, 0.50343, 0.50399, &
      0.50446/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction9U = &
    (/0.49235, 0.49301, 0.49375, 0.49437, &
      0.49485, 0.49531, 0.49569, 0.49594, 0.49612, &
      0.49623, 0.49631, 0.49637, 0.49641, 0.49645, &
      0.49651, 0.49658, 0.49669, 0.49684, 0.49703, &
      0.49728, 0.49754, 0.49779, 0.49807, 0.49837, &
      0.49866/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction15U = &
    (/0.04250, 0.04200, 0.04150, 0.04110, &
      0.04080, 0.04040, 0.04010, 0.03990, 0.03980, &
      0.03970, 0.03960, 0.03960, 0.03950, 0.03950, &
      0.03940, 0.03940, 0.03930, 0.03910, 0.03890, &
      0.03860, 0.03830, 0.03790, 0.03750, 0.03700, &
      0.03650/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction16U = &
    (/0.19600, 0.19130, 0.18710, 0.18330, &
      0.18040, 0.17730, 0.17480, 0.17290, 0.17160, &
      0.17060, 0.17000, 0.16950, 0.16920, 0.16890, &
      0.16840, 0.16780, 0.16690, 0.16560, 0.16370, &
      0.16110, 0.15830, 0.15510, 0.15160, 0.14680, &
      0.14280/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction17U = &
    (/0.57500, 0.58080, 0.58650, 0.59150, &
      0.59550, 0.59970, 0.60320, 0.60560, 0.60730, &
      0.60850, 0.60940, 0.61000, 0.61050, 0.61090, &
      0.61150, 0.61240, 0.61360, 0.61540, 0.61780, &
      0.62140, 0.62510, 0.62910, 0.63410, 0.63990, &
      0.64600/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction18U = &
    (/0.05350, 0.05320, 0.05300, 0.05270, &
      0.05260, 0.05240, 0.05220, 0.05210, 0.05200, &
      0.05200, 0.05200, 0.05190, 0.05190, 0.05190, &
      0.05190, 0.05180, 0.05180, 0.05170, 0.05160, &
      0.05140, 0.05130, 0.05110, 0.05090, 0.05060, &
      0.05030/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction19U = &
    (/0.21500, 0.21020, 0.20590, 0.20200, &
      0.19890, 0.19600, 0.19330, 0.19140, 0.19010, &
      0.18910, 0.18840, 0.18800, 0.18760, 0.18730, &
      0.18680, 0.18620, 0.18530, 0.18400, 0.18190, &
      0.17940, 0.17630, 0.17330, 0.16960, 0.16490, &
      0.16040/)
    real(r8), dimension(25), parameter :: vlimbSidebandFraction20U = &
    (/0.58130, 0.58690, 0.59350, 0.59800, &
      0.60200, 0.60630, 0.60960, 0.61210, 0.61380, &
      0.61500, 0.61590, 0.61650, 0.61700, 0.61740, &
      0.61800, 0.61880, 0.62010, 0.62190, 0.62430, &
      0.62800, 0.63150, 0.63550, 0.64070, 0.64620, &
      0.65240/)
    real(r8), dimension(129), parameter :: vlimbSidebandFraction23U = &
    (/0.46025, 0.46024, 0.46024, 0.46024, &
      0.46023, 0.46023, 0.46023, 0.46022, 0.46022, &
      0.46022, 0.46021, 0.46021, 0.46021, 0.46020, &
      0.46020, 0.46020, 0.46019, 0.46019, 0.46019, &
      0.46018, 0.46018, 0.46018, 0.46017, 0.46017, &
      0.46017, 0.46016, 0.46016, 0.46016, 0.46015, &
      0.46015, 0.46015, 0.46014, 0.46014, 0.46014, &
      0.46013, 0.46013, 0.46013, 0.46012, 0.46012, &
      0.46012, 0.46011, 0.46011, 0.46011, 0.46010, &
      0.46010, 0.46010, 0.46009, 0.46009, 0.46009, &
      0.46008, 0.46008, 0.46008, 0.46007, 0.46007, &
      0.46007, 0.46006, 0.46006, 0.46006, 0.46005, &
      0.46005, 0.46005, 0.46004, 0.46004, 0.46004, &
      0.46003, 0.46003, 0.46003, 0.46002, 0.46002, &
      0.46002, 0.46001, 0.46001, 0.46001, 0.46000, &
      0.46000, 0.46000, 0.45999, 0.45999, 0.45999, &
      0.45998, 0.45998, 0.45998, 0.45997, 0.45997, &
      0.45997, 0.45996, 0.45996, 0.45996, 0.45995, &
      0.45995, 0.45995, 0.45994, 0.45994, 0.45994, &
      0.45994, 0.45994, 0.45994, 0.45993, 0.45993, &
      0.45993, 0.45992, 0.45992, 0.45991, 0.45991, &
      0.45991, 0.45990, 0.45990, 0.45990, 0.45989, &
      0.45989, 0.45989, 0.45988, 0.45988, 0.45988, &
      0.45987, 0.45987, 0.45987, 0.45986, 0.45986, &
      0.45986, 0.45985, 0.45985, 0.45985, 0.45984, &
      0.45984, 0.45984, 0.45983, 0.45983, 0.45983/)
    real(r8), dimension(129), parameter :: vlimbSidebandFraction24U = &
    (/0.48014, 0.48014, 0.48014, 0.48014, &
      0.48014, 0.48014, 0.48014, 0.48014, 0.48014, &
      0.48014, 0.48014, 0.48014, 0.48014, 0.48014, &
      0.48014, 0.48014, 0.48014, 0.48015, 0.48015, &
      0.48015, 0.48015, 0.48015, 0.48015, 0.48015, &
      0.48015, 0.48015, 0.48015, 0.48015, 0.48015, &
      0.48015, 0.48015, 0.48015, 0.48015, 0.48015, &
      0.48015, 0.48015, 0.48015, 0.48015, 0.48015, &
      0.48015, 0.48016, 0.48016, 0.48016, 0.48016, &
      0.48016, 0.48016, 0.48016, 0.48016, 0.48016, &
      0.48016, 0.48016, 0.48016, 0.48016, 0.48016, &
      0.48016, 0.48016, 0.48016, 0.48016, 0.48016, &
      0.48016, 0.48016, 0.48016, 0.48016, 0.48017, &
      0.48017, 0.48017, 0.48017, 0.48017, 0.48017, &
      0.48017, 0.48017, 0.48017, 0.48017, 0.48017, &
      0.48017, 0.48017, 0.48017, 0.48017, 0.48017, &
      0.48017, 0.48017, 0.48017, 0.48017, 0.48017, &
      0.48017, 0.48017, 0.48018, 0.48018, 0.48018, &
      0.48018, 0.48018, 0.48018, 0.48018, 0.48018, &
      0.48018, 0.48018, 0.48018, 0.48018, 0.48018, &
      0.48018, 0.48018, 0.48018, 0.48018, 0.48018, &
      0.48018, 0.48018, 0.48018, 0.48018, 0.48018, &
      0.48019, 0.48019, 0.48019, 0.48019, 0.48019, &
      0.48019, 0.48019, 0.48019, 0.48019, 0.48019, &
      0.48019, 0.48019, 0.48019, 0.48019, 0.48019, &
      0.48019, 0.48019, 0.48019, 0.48019, 0.48019/)
    real(r8), dimension(129), parameter :: vlimbSidebandFraction25U = &
    (/0.49647, 0.49647, 0.49647, 0.49647, &
      0.49647, 0.49647, 0.49647, 0.49647, 0.49647, &
      0.49647, 0.49646, 0.49646, 0.49646, 0.49646, &
      0.49646, 0.49646, 0.49646, 0.49646, 0.49646, &
      0.49646, 0.49646, 0.49646, 0.49646, 0.49646, &
      0.49645, 0.49645, 0.49645, 0.49645, 0.49645, &
      0.49645, 0.49645, 0.49645, 0.49645, 0.49645, &
      0.49645, 0.49645, 0.49645, 0.49644, 0.49644, &
      0.49644, 0.49644, 0.49644, 0.49644, 0.49644, &
      0.49644, 0.49644, 0.49644, 0.49644, 0.49644, &
      0.49644, 0.49643, 0.49643, 0.49643, 0.49643, &
      0.49643, 0.49643, 0.49643, 0.49643, 0.49643, &
      0.49643, 0.49643, 0.49643, 0.49643, 0.49642, &
      0.49642, 0.49642, 0.49642, 0.49642, 0.49642, &
      0.49642, 0.49642, 0.49642, 0.49642, 0.49642, &
      0.49642, 0.49642, 0.49641, 0.49641, 0.49641, &
      0.49641, 0.49641, 0.49641, 0.49641, 0.49641, &
      0.49641, 0.49641, 0.49641, 0.49641, 0.49641, &
      0.49640, 0.49640, 0.49640, 0.49640, 0.49640, &
      0.49640, 0.49640, 0.49640, 0.49640, 0.49640, &
      0.49640, 0.49640, 0.49640, 0.49640, 0.49639, &
      0.49639, 0.49639, 0.49639, 0.49639, 0.49639, &
      0.49639, 0.49639, 0.49639, 0.49639, 0.49639, &
      0.49639, 0.49639, 0.49638, 0.49638, 0.49638, &
      0.49638, 0.49638, 0.49638, 0.49638, 0.49638, &
      0.49638, 0.49638, 0.49638, 0.49638, 0.49638/)
    real(r8), dimension(11), parameter :: vlimbSidebandFraction27U = &
    (/0.49004, 0.48991, 0.48981, 0.48973, &
      0.48966, 0.48962, 0.48957, 0.48949, 0.48939, &
      0.48923, 0.48897/)
    real(r8), dimension(4), parameter :: vlimbSidebandFraction33U = &
    (/0.48061, 0.48947, 0.50383, 0.50247/)

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
       use CFM_Vector_m, only: Vector_T, VectorValue_T
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
