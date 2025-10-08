import math
from typing import List, Optional
from rich.table import Table
from rich.console import Console
from rich.prompt import Prompt

def prime_factors(n: int) -> List[int]:
    """Compute the prime factors of an integer n."""
    factors: List[int] = []
    i = 2
    while i * i <= n:
        while n % i == 0:
            factors.append(i)
            n //= i
        i += 1
    if n > 1 or not factors:
        factors.append(n)
    return factors

def next_valid_lat(nlat: int, nsph: int, lat_mult: int) -> int:
    """Find the smallest multiple of nlat satisfying the anti-aliasing condition."""
    while 2 * nlat < 3 * (nsph - 1) + 1:
        nlat += lat_mult
    return nlat

def next_valid_lon(nlon: int, nfou: int, lon_maxprime: int) -> int:
    """Find the smallest nlon meeting FFT factorization and anti-aliasing constraints."""
    while nlon < 3 * nfou + 1 or max(prime_factors(nlon)) > lon_maxprime:
        nlon += 1
    return nlon


def get_grid_for_truncation(
        nfou: int,
        lat_mult: int = 4,
        lon_maxprime: int = 2,
    ) -> dict:
    """Return grid parameters for a specific truncation number."""
    nlat, nlon = lat_mult, 1
    nsph = nfou + 1
    nlat = next_valid_lat(nlat, nsph, lat_mult)
    nlon = next_valid_lon(nlon, nfou, lon_maxprime)
    size = nlat * nlon
    min_size = math.ceil((3 * (nsph - 1) + 1) / 2) * (3 * nfou + 1)
    over = size - min_size
    percent_over = "  - " if over == 0 else f"{100 * over / size:4.1f}"
    truncation_summary = {
        "nfou": nfou,
        "nsph": nsph,
        "nlat": nlat,
        "nlon": nlon,
        "over%": percent_over,
        "lon_factors": prime_factors(nlon),
    }
    return truncation_summary


def generate_grids_table(
        nfou_max: Optional[int] = None,
        lat_mult: int = 4,
        lon_maxprime: int = 2,
        latlon_maxprod: int = 1024 * 1024,
    ) -> Table:

    """Generate a rich.Table for spectral transform grids."""
    table = Table(title="Spectral Transform Grids")
    table.add_column("T")
    table.add_column("fou")
    table.add_column("sph")
    table.add_column("lat")
    table.add_column("lon")
    table.add_column("over%")
    table.add_column("lon factorization")

    nlat, nlon, nfou = lat_mult, 1, 1
    console = Console()

    while True:
        if nfou_max and nfou > nfou_max:
            break
        nsph = nfou + 1
        nlat = next_valid_lat(nlat, nsph, lat_mult)
        nlon = next_valid_lon(nlon, nfou, lon_maxprime)
        size = nlat * nlon
        if size > latlon_maxprod:
            exceeded = True
            break
        min_size = math.ceil((3 * (nsph - 1) + 1) / 2) * (3 * nfou + 1)
        over = size - min_size
        percent_over = "  - " if over == 0 else f"{100 * over / size:4.1f}"
        table.add_row(
            f"[bold cyan]T{nfou}[/bold cyan]",
            str(nfou),
            str(nsph),
            str(nlat),
            str(nlon),
            percent_over,
            ", ".join(map(str, prime_factors(nlon))),
        )
        nfou += 1
    return table, exceeded, size, latlon_maxprod



def main():
    console = Console()
    console.print(
        "[bold]Summary of Spectral Grid Sizes in Isca[/bold]\n"
        "This script will show you the minimum lat-lon grid shapes used in a\n"
        "given triangular truncation that satisfies both anti-aliasing and\n"
        "any (lat) multiple or (lon) factorisation constraints.\n\n"


        "[bold]Spectral Grid Summary[/bold]\n"
        "This table shows the spectral transform grids for some typical\n"
        "resolutions (e.g. T42, T85), with columns defined by:\n"
        "- [teal]T[/teal]: Truncation number\n"
        "- [teal]fou[/teal]: Fourier wavenumber\n"
        "- [teal]sph[/teal]: Spherical harmonic degree\n"
        "- [teal]lat[/teal]: Number of latitude grid points\n"
        "- [teal]lon[/teal]: Number of longitude grid points\n"
        "- [teal]over%[/teal]: Percentage over minimum required grid size required\n"
        "       to only satisfy anti-aliasing constraint\n"
        "- [teal]lon factorization[/teal]: Prime factors of longitude resolution.\n"
    )

    key_resolutions = [21, 42, 85, 170]
    table = Table(title="Common spectral transform grids (see below for more)")
    table.add_column("T")
    table.add_column("fou")
    table.add_column("sph")
    table.add_column("lat")
    table.add_column("lon")
    table.add_column("over%")
    table.add_column("lon factorization")

    for nfou in key_resolutions:
        grid = get_grid_for_truncation(nfou)
        table.add_row(
            f"[bold cyan]T{nfou}[/bold cyan]",
            str(grid["nfou"]),
            str(grid["nsph"]),
            str(grid["nlat"]),
            str(grid["nlon"]),
            grid["over%"],
            ", ".join(map(str, grid["lon_factors"])),
        )

    console.print(table)

    console.print(
        "\n[bold]Options:[/bold]\n"
        "a) Display full table of resolutions for all truncation numbers (T1 to Tx)\n"
        "b) Print grid parameters for a single specified truncation number (e.g. T21, T42)\n"
        "q) [yellow]Quit[/yellow]"
    )

    choice = Prompt.ask(
        "\nEnter your choice",
        choices=["a", "b", "q"],
        default="q",
    )

    if choice == "a":
        table, exceeded, size, latlon_maxprod = generate_grids_table()
        console.print(
            "\n[bold]Full Spectral Transform Grids Table[/bold]\n"
            "This table shows all available resolutions for spectral transform grids.\n"
        )
        console.print(table)
        if exceeded:
            console.print(f"[bold green]Next grid size ({size}) exceeds specified maximum {latlon_maxprod}[/bold green]")
    elif choice == "b":
        nfou = int(Prompt.ask("\nEnter truncation number (e.g., 21, 42, 85, 170)"))
        grid = get_grid_for_truncation(nfou)
        console.print(
            f"\n[bold]Grid parameters for T{nfou}:[/bold]\n"
            f"- Truncation number (T): {grid['nfou']}\n"
            f"- Spherical harmonic degree (sph): {grid['nsph']}\n"
            f"- Number of latitude grid points (lat): {grid['nlat']}\n"
            f"- Number of longitude grid points (lon): {grid['nlon']}\n"
            f"- Percentage over minimum grid size (over%): {grid['over%']}\n"
            f"- Longitude prime factors: {', '.join(map(str, grid['lon_factors']))}\n"
        )

if __name__ == "__main__":
    main()

