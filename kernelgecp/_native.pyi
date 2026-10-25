"""Optional private 128-bit reference backend; absent from universal wheels."""

NATIVE_PROTOCOL_VERSION: int
WORKING_PRECISION_BITS: int

def build_info() -> dict[str, str | int | bool]: ...

class ReferenceCardinal128:
    def __init__(
        self,
        precision_bits: int,
        omega_nodes: list[float],
        ordering: list[int],
        lower: list[list[str]],
        upper: list[list[str]],
    ) -> None: ...
    @property
    def rank(self) -> int: ...
    def design_tau(self, normalized_times: list[float]) -> bytes: ...
