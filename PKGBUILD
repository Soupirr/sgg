# Maintainer: Soupirr
pkgname=sgg
pkgver=0.1.0
pkgrel=1
pkgdesc="Identify pathogen genotypes and predict pathogenicity from nucleotide sequences"
arch=('any')
url="https://github.com/Soupirr/sgg"
license=('MIT')
depends=(
    'python'
    'python-pandas'
    'python-plotly'
    'python-biopython'
    'mafft'
    'fasttree'
    'iqtree'
)
makedepends=(
    'python-build'
    'python-installer'
    'python-wheel'
    'python-setuptools'
)
source=("$pkgname-$pkgver.tar.gz::https://github.com/Soupirr/sgg/archive/refs/tags/v$pkgver.tar.gz")
sha256sums=('b751e59bdabbce326c947ab07b928b8b8e7c272f75bc05a99bed6a4377d61c76')

build() {
    cd "$pkgname-$pkgver"
    /usr/bin/python -m build --wheel --no-isolation
}

package() {
    cd "$pkgname-$pkgver"
    /usr/bin/python -m installer --destdir="$pkgdir" dist/*.whl
    install -Dm644 LICENSE "$pkgdir/usr/share/licenses/$pkgname/LICENSE"
}
