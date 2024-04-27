from kiraan_waktu_solat import Takwim


class Pulau_Pinang(Takwim):
    def __init__(self):
        super().__init__()

    def takwim_solat_tahunan_multipoint_Penang(self, directory):
        super().takwim_solat_tahunan_multipoint(
            directory=directory,
            A=(5.282574, 100.177892, 40), B=(5.284480, 100.17720, 40),
            C=(5.142472, 100.502861,  40), D=(5.551423, 100.528866, 40),
            E=(5.472839, 100.180822, 40)
        )


class Perlis(Takwim):
    def __init__(self):
        super().__init__()

    def takwim_solat_tahunan_multipoint_Perlis(self, directory):
        super().takwim_solat_tahunan_multipoint(
            directory=directory,
            A=(6.422291, 100.120985, 10), B=(6.576966, 100.157826, 40),
            C=(6.692883, 100.172916, 40), D=(6.664744, 100.326753, 40),
            E=(6.531021, 100.370315, 40), F=(6.454471, 100.363438, 40),
            G=(6.254233, 100.193477, 10), H=(6.397214, 100.126859, 10)
        )


class Selangor(Takwim):
    def __init__(self):
        super().__init__()

    def takwim_solat_tahunan_multipoint_zon_1(self, directory):
        super().takwim_solat_tahunan_multipoint(
            directory=directory, KgGedangsA=(3.733333, 101.383333, 0),
            TgRhu=(2.637778, 101.616667, 0), PknBrogA=(2.94, 101.911389, 0),
            Titik1=(2.88079167, 101.87649722, 0),
            Titik3=(2.96680833, 101.58334167, 0),
            Titik4=(3.25626111, 101.85685833, 0),
            Titik6=(3.77920278, 101.35279444, 0)
            )

    def takwim_solat_tahunan_multipoint_kausar(self, directory):
        super().takwim_solat_tahunan_multipoint(
            directory=directory, Titik1=(2.88079167, 101.87649722, 0),
            Titik2=(2.61532500, 101.69205556, 0),
            Titik3=(2.96680833, 101.58334167, 0),
            Titik4=(3.25626111, 101.85685833, 0),
            Titik5=(3.61107500, 101.78599722, 0),
            Titik6=(3.77920278, 101.35279444, 0),
            Titik7=(2.94000000, 101.91138889, 0)
        )

    def takwim_solat_tahunan_multipoint_zon_2(self, directory):
        super().takwim_solat_tahunan_multipoint(
            directory=directory, BalaiCerapSelangor=(3.933333, 100.7, 10),
            PekanBrogA=(3.2025, 101.488611, 0)
        )

    def takwim_solat_tahunan_multipoint_zon_3(self, directory):
        super().takwim_solat_tahunan_multipoint(
            directory=directory, PulauKetam=(3.016667, 101.25, 0),
            TmnLangatMurni=(2.802778, 101.655, 40)
            )

    def takwim_solat_tahunan_multipoint_zon1_cadanganIzzat(self, directory):
        super().takwim_solat_tahunan_multipoint(
            directory=directory, titik1=(3.785131, 101.326643, 0),
            titik2=(3.598409, 101.349989, 0),
            titik3=(3.093907, 101.444230, 0),
            titik4=(2.642998, 101.594092, 0),
            titik5=(3.028763, 101.969515, 0),
            titik6=(3.136410, 101.959902, 0),
            titik7=(3.247474, 101.928317, 0),
            titik8=(3.607318, 101.808840, 0)
        )
