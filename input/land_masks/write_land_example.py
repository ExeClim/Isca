from isca import write_land

# create an Experiment object to handle the configuration of model parameters
# and output diagnostics

# Isca Continent Boundaries

write_land(resolution=42,land_mode='continents',continents=['all'],fname='Isca_continents_t42')
write_land(resolution=85,land_mode='continents',continents=['all'],fname='Isca_continents_t85')
write_land(resolution=170,land_mode='continents',continents=['all'],fname='Isca_continents_t170')

# 2x2 Island Archipelago

write_land(resolution=42,land_mode='square',boundaries=[[3,7,173,177],[3,7,183,187],\
                                                        [-7,-3,173,177],[-7,-3,183,187]],fname='2x2Islands_t42')
write_land(resolution=85,land_mode='square',boundaries=[[3,7,173,177],[3,7,183,187],\
                                                        [-7,-3,173,177],[-7,-3,183,187]],fname='2x2Islands_t85')
write_land(resolution=170,land_mode='square',boundaries=[[3,7,173,177],[3,7,183,187],\
                                                        [-7,-3,173,177],[-7,-3,183,187]],fname='2x2Islands_t170')

# 3x3 Island Archipelago

write_land(resolution=42,land_mode='square',boundaries=[[3,7,173,177],[3,7,183,187],[-7,-3,173,177],\
                                                        [-7,-3,183,187],[-2,2,178,182]],fname='3x3Islands_t42')
write_land(resolution=85,land_mode='square',boundaries=[[3,7,173,177],[3,7,183,187],[-7,-3,173,177],\
                                                        [-7,-3,183,187],[-2,2,178,182]],fname='3x3Islands_t85')
write_land(resolution=170,land_mode='square',boundaries=[[3,7,173,177],[3,7,183,187],[-7,-3,173,177],\
                                                        [-7,-3,183,187],[-2,2,178,182]],fname='3x3Islands_t170')

# 5x5 Island Archipelago

write_land(resolution=42,land_mode='square',boundaries=[[8,12,168,172],[8,12,178,182],[8,12,188,192],\
                                                        [3,7,173,177],[3,7,183,187],[-7,-3,173,177],[-7,-3,183,187],\
                                                        [-2,2,168,172],[-2,2,178,182],[-2,2,188,192],\
                                                        [-12,-8,168,172],[-12,-8,178,182],[-12,-8,188,192]],fname='5x5Islands_t42')
write_land(resolution=85,land_mode='square',boundaries=[[8,12,168,172],[8,12,178,182],[8,12,188,192],\
                                                        [3,7,173,177],[3,7,183,187],[-7,-3,173,177],[-7,-3,183,187],\
                                                        [-2,2,168,172],[-2,2,178,182],[-2,2,188,192],\
                                                        [-12,-8,168,172],[-12,-8,178,182],[-12,-8,188,192]],fname='5x5Islands_t85')
write_land(resolution=170,land_mode='square',boundaries=[[8,12,168,172],[8,12,178,182],[8,12,188,192],\
                                                        [3,7,173,177],[3,7,183,187],[-7,-3,173,177],[-7,-3,183,187],\
                                                        [-2,2,168,172],[-2,2,178,182],[-2,2,188,192],\
                                                        [-12,-8,168,172],[-12,-8,178,182],[-12,-8,188,192]],fname='5x5Islands_t170')

# Two Continents (Uniform Size / Symmetrical Position)

write_land(resolution=42,land_mode='square',boundaries=[[15,45,160,200],[-45,-15,160,200]],fname='2cont_unisym_t42')
write_land(resolution=85,land_mode='square',boundaries=[[15,45,160,200],[-45,-15,160,200]],fname='2cont_unisym_t85')
write_land(resolution=170,land_mode='square',boundaries=[[15,45,160,200],[-45,-15,160,200]],fname='2cont_unisym_t170')

# Two Continents (Large + Small / Symmetrical Position)

write_land(resolution=42,land_mode='square',boundaries=[[15,75,135,225],[-45,-15,160,200]],fname='2cont_LSsym_t42')
write_land(resolution=85,land_mode='square',boundaries=[[15,75,135,225],[-45,-15,160,200]],fname='2cont_LSsym_t85')
write_land(resolution=170,land_mode='square',boundaries=[[15,75,135,225],[-45,-15,160,200]],fname='2cont_LSsym_t170')

# Two Continents (Large + Small / Asymmetrical Position)

write_land(resolution=42,land_mode='square',boundaries=[[15,75,135,225],[-45,-15,185,225]],fname='2cont_LSasym_t42')
write_land(resolution=85,land_mode='square',boundaries=[[15,75,135,225],[-45,-15,185,225]],fname='2cont_LSasym_t85')
write_land(resolution=170,land_mode='square',boundaries=[[15,75,135,225],[-45,-15,185,225]],fname='2cont_LSasym_t170')

# Two Continents (UniSym) + 2x2 Island Archipelago

write_land(resolution=42,land_mode='square',boundaries=[[3,7,173,177],[3,7,183,187],\
                                                        [-7,-3,173,177],[-7,-3,183,187],\
                                                        [15,45,160,200],[-45,-15,160,200]],fname='2contUniSym_2x2Islands_t42')
write_land(resolution=85,land_mode='square',boundaries=[[3,7,173,177],[3,7,183,187],\
                                                        [-7,-3,173,177],[-7,-3,183,187],\
                                                        [15,45,160,200],[-45,-15,160,200]],fname='2contUniSym_2x2Islands_t85')
write_land(resolution=170,land_mode='square',boundaries=[[3,7,173,177],[3,7,183,187],\
                                                        [-7,-3,173,177],[-7,-3,183,187],\
                                                        [15,45,160,200],[-45,-15,160,200]],fname='2contUniSym_2x2Islands_t170')

# Two Continents (UniSym) + 3x3 Island Archipelago

write_land(resolution=42,land_mode='square',boundaries=[[3,7,173,177],[3,7,183,187],[-7,-3,173,177],\
                                                        [-7,-3,183,187],[-2,2,178,182],\
                                                        [15,45,160,200],[-45,-15,160,200]],fname='2contUniSym_3x3Islands_t42')
write_land(resolution=85,land_mode='square',boundaries=[[3,7,173,177],[3,7,183,187],[-7,-3,173,177],\
                                                        [-7,-3,183,187],[-2,2,178,182],\
                                                        [15,45,160,200],[-45,-15,160,200]],fname='2contUniSym_3x3Islands_t85')
write_land(resolution=170,land_mode='square',boundaries=[[3,7,173,177],[3,7,183,187],[-7,-3,173,177],\
                                                        [-7,-3,183,187],[-2,2,178,182],\
                                                        [15,45,160,200],[-45,-15,160,200]],fname='2contUniSym_3x3Islands_t170')

# Two Continents (UniSym) + 5x5 Island Archipelago

write_land(resolution=42,land_mode='square',boundaries=[[8,12,168,172],[8,12,178,182],[8,12,188,192],\
                                                        [3,7,173,177],[3,7,183,187],[-7,-3,173,177],[-7,-3,183,187],\
                                                        [-2,2,168,172],[-2,2,178,182],[-2,2,188,192],\
                                                        [-12,-8,168,172],[-12,-8,178,182],[-12,-8,188,192],\
                                                        [15,45,160,200],[-45,-15,160,200]],fname='2contUniSym_5x5Islands_t42')
write_land(resolution=85,land_mode='square',boundaries=[[8,12,168,172],[8,12,178,182],[8,12,188,192],\
                                                        [3,7,173,177],[3,7,183,187],[-7,-3,173,177],[-7,-3,183,187],\
                                                        [-2,2,168,172],[-2,2,178,182],[-2,2,188,192],\
                                                        [-12,-8,168,172],[-12,-8,178,182],[-12,-8,188,192],\
                                                        [15,45,160,200],[-45,-15,160,200]],fname='2contUniSym_5x5Islands_t85')
write_land(resolution=170,land_mode='square',boundaries=[[8,12,168,172],[8,12,178,182],[8,12,188,192],\
                                                        [3,7,173,177],[3,7,183,187],[-7,-3,173,177],[-7,-3,183,187],\
                                                        [-2,2,168,172],[-2,2,178,182],[-2,2,188,192],\
                                                        [-12,-8,168,172],[-12,-8,178,182],[-12,-8,188,192],\
                                                        [15,45,160,200],[-45,-15,160,200]],fname='2contUniSym_5x5Islands_t170')

