from math import *
from errors import *


class _rod:
    def __init__(self):
        self.temperature = 100
        self.old_flux = 0
        self.flux = 0
        self.fuel = 1
        self.delta_flux = 0


class _reactor:
    def __init__(self, depth, size):
        self.apr = 0
        self.iodine = 0
        self.xenon = 0
        self.water_density = 1000
        self.max_neutrons = depth * size**2 * 10**8
        center_offset = (size - 1) / 2

        self.rod_values = [
            [
                (
                    1
                    if sqrt((y - center_offset) ** 2 + (x - center_offset) ** 2)
                    < size / 2
                    else None
                )
                for y in range(size)
            ]
            for x in range(size)
        ]

        self.cells = [
            [
                [
                    _rod() if self.rod_values[x][y] is not None else None
                    for y in range(size)
                ]
                for x in range(size)
            ]
            for z in range(depth)
        ]

    def update(self):
        for depth in self.cells:
            for row in depth:
                for cell in row:
                    self.apr += cell.flux
        self.apr /= self.max_neutrons

        iodine_scaling = 0.00025
        self.iodine += self.apr * iodine_scaling
        iodine_decay = self.iodine * iodine_scaling
        self.iodine -= iodine_decay

        xenon_timefactor = 1
        self.xenon += iodine_decay
        self.xenon -= 2 / 3 * self.xenon * iodine_scaling * xenon_timefactor
        self.xenon -= self.xenon * self.apr * iodine_scaling

        spreading_factors = [
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [-1, 0, 0],
            [-1, -1, 0],
            [-1, 0, 0],
            [-1, 1, 0],
            [1, -1, 0],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
            [-1, 0, 1],
            [-1, -1, 1],
            [-1, 0, 1],
            [-1, 1, 1],
            [1, -1, 1],
            [1, 0, -1],
            [1, 1, -1],
            [0, 1, -1],
            [-1, 0, -1],
            [-1, -1, -1],
            [-1, 0, -1],
            [-1, 1, -1],
            [1, -1, -1],
            [0, 0, 1],
            [0, 0, -1],
        ]

        for depth in self.cells:
            for row in depth:
                for cell in row:
                    cell.delta_flux = 0

        for depth, section in enumerate(self.cells):
            for x, row in enumerate(section):
                for y, cell in enumerate(row):
                    if cell is None:
                        return

                    flux_exp = 2 * (0.2 + cell.rod_percentage * 0.8)
                    flux_exp *= self.water_density / 1000
                    flux_exp *= cell.fuel
                    flux_exp *= 1 - self.xenon / 4.5
                    flux_exp *= 1 - 0.3 * self.apr
                    flux_exp /= 2  # neutron spread factor

                    cell.flux = flux_exp * (10 + cell.flux)

                    for factor in spreading_factors:
                        try:
                            other_cell = self.cells[depth + factor[0]][x + factor[1]][
                                y + factor[2]
                            ]
                        except IndexError:
                            continue

                        other_cell.delta_flux += (
                            cell.flux / 26 * (1 - cell.rod_percentage)
                        )

        for depth in self.cells:
            for row in depth:
                for cell in row:
                    cell.flux += cell.delta_flux

    def withdraw(self, locations, amount):
        height = len(self.cells)

        for location in locations:
            self.rod_values[location[0]][location[1]] -= amount
            remaining_rod = self.rod_values[location[0]][location[1]]
            for depth, section in enumerate(self.cells):
                if not remaining_rod:
                    break

                remaining_rod -= 1 / height
                if remaining_rod >= 1 / height:
                    new_percentage = 1
                else:
                    new_percentage = remaining_rod * height
                section[location[0]][location[1]].rod_percentage = new_percentage

    def insert(self, locations, amount):
        self.withdraw(locations, -amount)
