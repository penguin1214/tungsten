/*
 * Use RIGHT-HAND coordinate system.
 * All primitives normalize to ([0,1], [0.1], [0,1])
 * **ATTENTION** NOT CONSISTENT WITH PLAINGRID NOW!
 * SPP: sample per pixel
 */

/*
 * int_medium vs. ext_medium???
 */
#include "Version.hpp"
#include "Shared.hpp"
#include "io/MeshIO.hpp"

using namespace Tungsten;

int main(int argc, const char *argv[])
{
	CliParser parser("tungsten", "[options] scene1 [scene2 [scene3...]]");

	StandaloneRenderer renderer(parser, std::cout);

	parser.parse(argc, argv);

	if (parser.isPresent(OPT_VERSION)) {
		std::cout << "tungsten, version " << VERSION_STRING << std::endl;
		std::exit(0);
	}

	/// TODO
	renderer.setup();

	/// load scenes
	/// render
	while (renderer.renderScene());

	return 0;
}
