FROM rust:1.69

COPY ./ /mergebams

WORKDIR /mergebams

RUN cargo build --release

RUN chmod a+x /mergebams/target/release/mergebams

ENV PATH "/mergebams/target/release:$PATH"

CMD ["/bin/bash"]
